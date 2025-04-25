---
layout: archive
title: ""   
permalink: /CS/
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

 

##  <span style="color:#8C8C8C"> Compressive sensing</span> 
---
<div style="text-align:justify;">
Compressive sensing (CS) problems aim to recover a sparse signal $\mathbf{x}^*\in\mathbb{R}^{n}$ from the following linear system,
</div>

\begin{equation}
\mathbf{b} = \mathbf{A}\mathbf{x} + \boldsymbol{\varepsilon} \tag{CS}
\end{equation} 

<div style="text-align:justify;">
where $\mathbf{A}\in\mathbb{R}^{m\times n}$ is the sensing matrix, $\mathbf{b}\in\mathbb{R}^{m}$ is the observation, and $\boldsymbol{\varepsilon}\in\mathbb{R}^{n}$ is the noise. To recover the signal, the folowing optimzation models are freguently explored:
</div>   

<p style="line-height: 2;"></p>
<div style="text-align:justify;">
◻️ Sparsity constrained CS
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n}}~ \frac{1}{2}\parallel\mathbf{A}\mathbf{x}-\mathbf{b} \parallel^2,~~~\textrm{s.t.}~ \parallel\mathbf{x} \parallel_0\leq s \tag{SCCS}
\end{equation}
</div> 
<div style="text-align:justify;">
◻️ $L_q, 0\leq q <1$ norm regularized CS
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n}}~ \frac{1}{2}\parallel\mathbf{A}\mathbf{x}-\mathbf{b} \parallel^2+\lambda \parallel\mathbf{x} \parallel_q^q \tag{LqRCS}
\end{equation}
</div>   
<div style="text-align:justify;">
◻️ Reweighted L1 norm regularized CS
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n}}~ \frac{1}{2}\parallel\mathbf{A}\mathbf{x}-\mathbf{b} \parallel^2+\lambda \parallel \mathbf{W} \mathbf{x} \parallel_1 \tag{RL1CS}
\end{equation} 
where $\parallel\mathbf{x}\parallel_q^q=\sum_i |x_i|^q$ denotes the $\ell_q$-norm, in particular, $\parallel\mathbf{x}\parallel_0:=\parallel\mathbf{x}\parallel_0^0$ denotes the $\ell_0$-norm, which counts the number of nonzero entries in $\mathbf{x}$,  $\lambda>0$ is the penalty parameter, and $\mathbf{W}$ is a diagonal matrix with positive diagonal entrices.
</div> 
---
<div style="text-align:justify;">
The package can be download here - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href="\files\CSpack.zip" target="_blank">CSpack</a>, which provides 6 solvers from the following papers, where <b style="font-size:16px;color:#777777">NHTP</b>, <b style="font-size:16px;color:#777777">GPSP</b>, and <b style="font-size:16px;color:#777777">IIHT</b> are designed to solve (SCCS), <b style="font-size:16px;color:#777777">PSNP</b> is designed to solve (LqRCS) with  $0\leq q < 1$,  <b style="font-size:16px;color:#777777">NL0R</b> is designed to solve (L0RCS) with  $q=0$, and <b style="font-size:16px;color:#777777">MIRL1</b> is designed to solve (RL1CS).
</div>  
> <b style="font-size:14px;color:#777777">NHTP</b> - <span style="font-size: 14px"> S Zhou, N Xiu, and H Qi, Global and quadratic convergence of Newton hard-thresholding pursuit, JMLR, 22:1-45, 2021. </span>
<br><b style="font-size:14px;color:#777777">GPSP</b> - <span style="font-size: 14px"> S Zhou, Gradient projection Newton pursuit for sparsity constrained optimization, ACHA, 61:75-100, 2022. </span>
<br><b style="font-size:14px;color:#777777">PSNP</b> - <span style="font-size: 14px"> S Zhou, X Xiu, Y Wang, and D Peng, Revisiting Lq ( 0 ≤ q < 1 ) norm regularized optimization, arXiv:2306.14394, 2023. </span>
<br><b style="font-size:14px;color:#777777">NL0R</b> - <span style="font-size: 14px"> S Zhou, L Pan, and N Xiu, Newton method for l0 regularized optimization, Numer Algorithms, 88:1541–1570, 2021. </span>
<br><b style="font-size:14px;color:#777777">IIHT</b> - <span style="font-size: 14px"> L Pan, S Zhou, N Xiu, and H Qi, A convergent iterative hard thresholding for nonnegative sparsity optimization, PJO, 13:325-353, 2017. </span>
<br><b style="font-size:14px;color:#777777">MIRL1</b> - <span style="font-size: 14px"> S Zhou, N Xiu, et. al, A Null-space-based weighted l1 minimization approach to compressed sensing, Inf inference, 5:76-102, 2016. </span>

---
<div style="text-align:justify;">
Below is a demonstration of how CSpack can be used to solve the problem. You need to input data $(\texttt{A},\texttt{At},\texttt{b},\texttt{n},\texttt{s})$  and select $\texttt{solver}$ from $\texttt{\{`NHTP',`GPNP',`IIHT',`PSNP',`NL0R',`MILR1'\}}$. 
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

solver  = {'NHTP', 'GPNP', 'IIHT', 'PSNP', 'NL0R', 'MIRL1'};
out     = CSpack(A,[],b,n,s,solver{1}); 

fprintf(' Objective of xopt:       %.2e\n', norm(A*xopt-b)^2/2);
fprintf(' Objective of out.sol:    %.2e\n', out.obj);
fprintf(' Sparsity of out.sol:     %2d\n', nnz(out.sol));
fprintf(' Computational time:      %.3fsec\n',out.time); 
```

<div style="text-align:justify;">
The inputs and outputs of CSpack are detailed below, where inputs $(\texttt{A},\texttt{At},\texttt{b},\texttt{n},\texttt{s},\texttt{solver})$ are required. If $\texttt{A}$ is a matrix,  $\texttt{At}$ can be $\texttt{A}'$ or $\texttt{[]}$. If $\texttt{A}$ is a function handle, then $\texttt{At}$ must be provided. If $\texttt{solver}$ is one of $\texttt{\{`PSNP',`NL0R',`MILR1'\}}$, then $\texttt{s}$ can be $\texttt{[]}$.  If $\texttt{solver}$ is one of $\texttt{\{`NHTP',`GPNP',`IIHT'\}}$, then $\texttt{s}$ must be provided. The parameters in $\texttt{pars}$ are optional, but setting certain ones can improve the solver's performance and the quality of the solution.
</div>

<p style="line-height: 1;"></p>

```ruby
function out = CSpack(A,At,b,n,s,solver,pars)
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
%
% Outputs:
%     out.sol:   The sparse solution x
%     out.sp:    Sparsity level of out.sol
%     out.time:  CPU time
%     out.iter:  Number of iterations
%     out.obj:   Objective function value at out.sol 
% -------------------------------------------------------------------------
```
