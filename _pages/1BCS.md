---
layout: archive
title: ""   
permalink: /1BCS/
author_profile: true
---

<style>
a:link {
  text-decoration: none;
}

a:visited {
  text-decoration: none;
}

a:hover {
  text-decoration: underline;
}

a:active {
  text-decoration: underline;
}
</style>

 

##  <span style="color:#8C8C8C"> One-bit compressive sensing</span> 
---
<div style="text-align:justify;">
One-bit compressive sensing (1BCS) problems aim to recover a signal $\mathbf{x}^*\in\mathbb{R}^{n}$ from the following system,
</div>

\begin{equation}
\mathbf{c} = \mathrm{Diag}(\mathbf{h}) \mathrm{sign}(\boldsymbol{\Phi}\mathbf{x} + \boldsymbol{\varepsilon}) \tag{1bCS}
\end{equation} 

<div style="text-align:justify;">
where $\boldsymbol{\Phi}\in\mathbb{R}^{m\times n}$ is the sensing matrix, both $\mathbf{c}$ and $\mathbf{h}\in\{-1,1\}^{m}$, $\boldsymbol{\varepsilon}\in\mathbb{R}^{n}$ is the noise, $\mathrm{Diag}(\mathbf{h})$ denotes the diagonal matrix with diagonal entries formed by $\mathbf{h}$, and $\mathrm{sign}(t)$ is the sign function of $t$ defined by $\mathrm{sign}(t)=1$ if $t>0$ and $\mathrm{sign}(t)=-1$ otherwise. Then $\mathrm{sign}(\mathbf{x})=(\mathrm{sign}(x_1),\ldots,\mathrm{sign}(x_n))^\top$. Note that multiplying $\mathrm{Diag}(\mathbf{h})$ means that the sign flips occurs when observed $\mathbf{c}$, making the problem harder. In this model, we assume at most $k$ signs are flipped, namely, $\mathbf{h}$ satisfies $\parallel\mathbf{h}-1\parallel_0\leq k$, where $k$ is a given scalar $\parallel\mathbf{x}\parallel_0$ denotes the $\ell_0$-norm, which counts the number of nonzero entries in $\mathbf{x}$. To recover the signal, the following double-sparsity constrained optimization (DSCO) model is solved, 
</div>      
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n},\mathbf{y}\in\mathbb{R}^{m}}~  \parallel\mathbf{A}\mathbf{x}+\mathbf{y} -\epsilon \parallel^2 + \eta \parallel \mathbf{x} \parallel^2,~~~\textrm{s.t.}~ \parallel\mathbf{x} \parallel_0\leq s,~ \parallel \mathbf{y}_+\parallel_0\leq k \tag{DSCO}
\end{equation}
<div style="text-align:justify;">
where $\mathbf{A} = \mathrm{Diag}(\mathbf{c}) \boldsymbol{\Phi}$, $\epsilon, \eta>0$ are two given scalar, and $\mathbf{y}_+=(\max\{0,y_1\},\ldots,\max\{0,y_m\})^\top$.
</div> 
  
---
<div style="text-align:justify;">
The solver can be download here - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href="https://github.com/ShenglongZhou/GPSP" target="_blank">GPSP</a>, which was developed from the following paper:
</div>  
>  <span style="font-size: 14px"> S Zhou, Z Luo, N Xiu, and G Li, Computing one-bit compressive sensing via double-sparsity constrained optimization, IEEE Tran Signal Process, 70:1593-1608, 2022. </span>
  
---
<div style="text-align:justify;">
The citation for GPSP is shown below. Here, inputs $(\texttt{A},\texttt{b},\texttt{n},\texttt{solver})$ are required. If $\texttt{A}$ is a function handle, then $\texttt{At}$ is required. If $\texttt{A}$ is a matrix,  $\texttt{At}$ can be $\texttt{A}'$ or $\texttt{[]}$. If $\texttt{solver}$ is one of $\texttt{\{`NHTP',`GPNP',`IIHT'\}}$, then $\texttt{s}$ is required. If $\texttt{solver}$ is one of $\texttt{\{`PSNP',`NL0R',`MILR1'\}}$, then $\texttt{s}$ can be $\texttt{[]}$.
</div>

<p style="line-height: 1;"></p>

```ruby
function out = CSsolver(A,At,b,n,s,solver,pars)
% -------------------------------------------------------------------------
% This solver solves compressive sensing (CS) in one of the following forms
%
% 1) The sparsity constrained compressive sensing (SCCS)
%
%         min_{x\in R^n} 0.5||Ax-b||^2  s.t. ||x||_0<=s
%
% 2) The L0 regularized compressive sensing (LqCS)
%
%         min_{x\in R^n} 0.5||Ax-b||^2 + lambda * ||x||_q^q,  0<=q<1 
%
% 3) The reweighted L1-regularized compressive sensing (RL1CS)
%
%         min_{x\in R^n} 0.5||Ax-b||^2 + lambda||Wx||_1
%
% where s << n is the given sparsity and lambda>0 
%       A\in\R{m by n} the measurement matrix
%       b\in\R{m by 1} the observation vector 
%       W\in\R{n by n} is a diagonal matrix to be updated iteratively
% -------------------------------------------------------------------------
% Inputs:
%   A  :     The measurement matrix, A\in\R{m by n}              (REQUIRED)
%   At :     The transpose of A and can be [] if A is a matrix   (REQUIRED)
%            But At is REQUIRED if A is a function handle 
%            i.e., A*x = A(x); A'*y = At(y); 
%   b:       The observation vector  b\in\R{m by 1}              (REQUIRED)
%   n:       Dimension of the solution x,                        (REQUIRED)
%   s:       The sparsity level, if unknown, put it as []        (REQUIRED)
%   solver:  A text string, can be one of                        (REQUIRED)
%            {'NHTP','GPNP','PSNP','NL0R','IIHT','MILR1'}
%
%           --------------------------------------------------------------------------------
%                    |  'NHTP'   |  'GPNP'   |  'PSNP'   |  'NL0R'   |  'IIHT'   |  'MIRL1'   
%           --------------------------------------------------------------------------------
%           Model    |   SCCS    |   SCCS    |   LqRCS   |   L0RCS   |   SCCS    |   RL1CS     
%           Method   | 2nd-order | 2nd-order | 2nd-order | 2nd-order | 1st-order | 1st-order  
%           Sparsity | required  | required  |  no need  |  no need  | required  |  no need
%           --------------------------------------------------------------------------------  
%
%   pars  : ----------------For all solvers -------------------------------
%           pars.x0     --  Starting point of x       (default, zeros(n,1))                     
%           pars.disp   --  =1, show results for each step      (default,1)
%                           =0, not show results for each step
%           pars.maxit  --  Maximum number of iterations     (default, 2e3) 
%           pars.tol    --  Tolerance of stopping criteria   (default,1e-6)
%           ----------------Particular for NHTP ---------------------------
%           pars.eta    --  A positive scalar for 'NHTP'       (default, 1)  
%                           Tuning pars.eta may improve solution quality.
%           ----------------Particular for PSNP ---------------------------
%           pars.q      --  Decide Lq norm                  (default,  0.5)  
%           pars.lambda --  An initial penalty parameter    (default,  0.1)
%           pars.obj    --  A predefined lower bound of f(x)(default,1e-20)
%           ----------------Particular for NL0R ---------------------------
%           pars.tau    --  A positive scalar for 'NL0R'    (default,    1)  
%           pars.lambda --  An initial penalty parameter    (default,  0.1)
%           pars.obj    --  A predefined lower bound of f(x)(default,1e-20)
%           ----------------Particular for IIHT ---------------------------
%           pars.neg    --  =0, Compute SCCS without x>=0       (default,0)
%                           =1, Compute SCCS with x>=0
% -------------------------------------------------------------------------
```

<div style="text-align:justify;">
Below is a demonstration of how CSpack can be used to solve the problem. You simply need to input the data $(\texttt{A},\texttt{b},\texttt{n},\texttt{s})$  and select $\texttt{solver}$ from $\texttt{\{`NHTP',`GPNP',`IIHT',`PSNP',`NL0R',`MILR1'\}}$. The parameters in $\texttt{pars}$ are optional, but setting certain ones can improve the solver's performance and the quality of the solution.
</div>

<p style="line-height: 1;"></p>

```ruby
clc; clear; close all; addpath(genpath(pwd));

n       = 10000;  
m       = ceil(0.25*n); 
s       = ceil(0.05*n); 

T       = randperm(n,s);  
xopt    = zeros(n,1);
xopt(T) = (0.1+rand(s,1)).*sign(randn(s,1));  
A       = randn(m,n)/sqrt(m);   
b       = A(:,T)*xopt(T)+0.00*randn(m,1);  

t       = 2; 
solver  = {'NHTP', 'GPNP', 'IIHT', 'PSNP', 'NL0R', 'MIRL1'};
out     = CSsolver(A,[],b,n,s,solver{t}); 

fprintf(' Objective of xopt:       %.2e\n', norm(A*xopt-b)^2/2);
fprintf(' Objective of out.sol:    %.2e\n', out.obj);
fprintf(' Sparsity of out.sol:     %2d\n', nnz(out.sol));
fprintf(' Computational time:      %.3fsec\n',out.time); 
```
