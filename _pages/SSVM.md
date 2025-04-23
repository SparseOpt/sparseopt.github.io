---
layout: archive
title: ""   
permalink: /SSVM/
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

 

##  <span style="color:#8C8C8C"> Sparse support vector mechine</span> 
---
<div style="text-align:justify;">
  The  soft-margin support vector mechine (sm-SVM) takes the form of 
</div>

\begin{equation}
\min_{(\mathbf{w};b)\in\mathbb{R}^{n+1}}~\frac{1}{2}\parallel \mathbf{w} \parallel^2 + C \sum_{i=1}^m\ell\left(1-y_i(b+  \mathbf{a}_i^\top\mathbf{w})\right) \tag{SVM}
\end{equation} 

<div style="text-align:justify;">
where $\mathbf{A}=(\mathbf{a}_1,\ldots,\mathbf{a}_m)\in\mathbb{R}^{n\times m}$ is the sample matrix, $\mathbf{y}=(y_1,\ldots,y_m)^\top\in\mathbb{R}^{m}$ is the label vector, $C>0$ is the penalty parameter, and $\ell$ is a loss function. One popular loss function is the hinge loss defined by  $\ell_{h}(t)=\max\{t,0\}.$ Two loss functions are considered below, resulting in two SVM models.
</div>      

<p style="line-height: 2;"></p>

◻️ Step function regularied SVM
\begin{equation}
\min_{(\mathbf{w};b)\in\mathbb{R}^{n+1}}~\frac{1}{2}\parallel \mathbf{w} \parallel^2 + C \sum_{i=1}^m\ell_{0/1}\left(1-y_i(b+  \mathbf{a}_i^\top\mathbf{w})\right) \tag{L01SVM}
\end{equation} 
<div style="text-align:justify;">
where $\ell_{0/1}$ is the step function (or 0/1 loss function) defined by $\ell_{0/1}(t)=1$ if $t>0$ and $\ell_{0/1}(t)=0$ otherwise. By letting $\mathbf{z}_+=(\max\{0,z_1\},$ $\ldots,$ $\max\{0,z_m\})^\top$ and $\parallel\mathbf{x}\parallel_0$ denote the $\ell_0$-norm, which counts the number of nonzero entries in $\mathbf{x}$, one can check that $ \sum_{i=1}^m\ell_{0/1}\left(1-y_i(b+  \mathbf{a}_i^\top\mathbf{w})\right) = \|  (1-\mathbf{A}\mathbf{w}-b\mathbf{y} )_+ \|_0$. 
</div>

<!--
◻️ $\ell_{cC}$ regularized  SVM
\begin{equation}
\min_{(\mathbf{w};b)\in\mathbb{R}^{n+1}}~\frac{1}{2} \parallel  \mathbf{w} \parallel^2 + \sum_{i=1}^m\ell_{cC}\left(1-y_i(b+  \mathbf{a}_i^\top\mathbf{w})\right) \tag{SFRSVM}
\end{equation} 
<div style="text-align:justify;">
where  $\ell_{cC}(t)=Ct^2/2$ if $t>0$ and $\ell_{cC}(t)=ct^2/2$ otherwise with $C>c>0$. The dual problem of (LcCSVM) is the following quadratic kernel-based SVM problem
</div>  

\begin{equation}
\min_{\boldsymbol{\alpha}\in\mathbb{R}^{m}}~\frac{1}{2} \boldsymbol{\alpha}^\top \mathbf{Q} \boldsymbol{\alpha} +\frac{1}{2}\sum_{i=1}^m h_{cC}(\alpha_i) -\mathbf{e}^\top\boldsymbol{\alpha}, ~~~~ \text{s.t.} ~~\mathbf{y}^\top\boldsymbol{\alpha}=0\tag{QKSVM}
\end{equation} 
<div style="text-align:justify;">
where $\mathbf{Q}=(Q_{ij})_{1\leq i,j\leq m}$ with $Q_{ij}=y_iy_j\mathbf{a}_i^\top\mathbf{a}_j$, $\mathbf{e}=(1,\ldots,1)^\top$, and $h_{cC}(t)=t^2/C$ if $t>0$ and $\ell_{cC}(t)=t^2/c$.
</div>  
-->

◻️ Sparsity constrained quadratic kernel-based SVM 
\begin{equation}
\min_{\boldsymbol{\alpha}\in\mathbb{R}^{m}}~\frac{1}{2} \boldsymbol{\alpha}^\top \mathbf{Q} \boldsymbol{\alpha} +\frac{1}{2}\sum_{i=1}^m h_{cC}(\alpha_i) -\mathbf{e}^\top\boldsymbol{\alpha}, ~~~~ \text{s.t.} ~~\mathbf{y}^\top\boldsymbol{\alpha}=0,~\parallel  \boldsymbol{\alpha} \parallel_0\leq s \tag{SCSVM}
\end{equation} 
<div style="text-align:justify;">
where $\mathbf{Q}\in\mathbb{R}^{m\times m}$ with each entry $Q_{ij}=y_iy_j\mathbf{a}_i^\top\mathbf{a}_j$, $\mathbf{e}=(1,\ldots,1)^\top$, and $h_{cC}(t)=t^2/C$ if $t>0$ and $\ell_{cC}(t)=t^2/c$,  $C>c>0$, and $s\ll m$. In fact, model (SCSVM) without sparsity constraint $\parallel  \boldsymbol{\alpha} \parallel_0\leq s$ is the dual problem of model (SVM) with $\ell=\ell_{cC}$, where  $\ell_{cC}(t)=t^2/2$ if $t>0$ and $\ell_{cC}(t)=(c/C)t^2/2$ otherwise. 
</div>  

> <div style="text-align:justify;"> According to the Representer Theorem,  otpimal soultion $ \mathbf{w}^* $ to (SVM) and otpimal soultion $\boldsymbol{\alpha}^* $ to the dual SVM satisfy $ \mathbf{w}^* = \sum_{i=1}^m \alpha_i^* y_i \mathbf{a}_i $. The training vectors $ \mathbf{a}_i $ corresponding to nonzero $ \alpha_i^* $ are known as support vectors. Therefore, both model (SFRSVM) and model (SCSVM) enable the reduction of support vectors. </div> 

---
<div style="text-align:justify;">
The package can be download here - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href=" " target="_blank">SSVMpack</a>, which provides 3 solvers from the following 3 papers, where <b style="font-size:16px;color:#777777">L01ADMM</b> and <b style="font-size:16px;color:#777777">NM01</b>, and <b style="font-size:16px;color:#777777">IIHT</b> are designed to solve (SFRSVM), and <b style="font-size:16px;color:#777777">NSSVM</b> is designed to solve (SCSNM).
</div>  

> <div style="text-align:justify;"> <b style="font-size:14px;color:#777777">NM01</b>-<span style="font-size: 14px"> S Zhou, L Pan, N Xiu, and H Qi, Quadratic convergence of smoothing Newton's method for 0/1 loss optimization, SIAM J Optim, 31:3184–3211, 2021. </span> </div>
> <div style="text-align:justify;">  <b style="font-size:14px;color:#777777">NSSVM</b>-<span style="font-size: 14px"> S Zhou, Sparse SVM for sufficient data reduction, IEEE Trans Pattern Anal  Mach  Intell, 44:5560-5571, 2022. </span> </div>
> <div style="text-align:justify;"> <b style="font-size:14px;color:#777777">L01ADMM</b>-<span style="font-size: 14px"> H Wang, Y Shao, S Zhou, C Zhang, and N Xiu, Support vector machine classifier via L0/1 soft-margin loss, IEEE Trans Pattern Anal  Mach  Intell, 44:7253-7265, 2022. </span> </div>

---
<div style="text-align:justify;">
The citation for CSpack is shown below. Here, inputs $(\texttt{A},\texttt{b},\texttt{n},\texttt{solver})$ are required. If $\texttt{A}$ is a function handle, then $\texttt{At}$ is required. If $\texttt{A}$ is a matrix,  $\texttt{At}$ can be $\texttt{A}'$ or $\texttt{[]}$. If $\texttt{solver}$ is one of $\texttt{\{`NHTP',`GPNP',`IIHT'\}}$, then $\texttt{s}$ is required. If $\texttt{solver}$ is one of $\texttt{\{`PSNP',`NL0R',`MILR1'\}}$, then $\texttt{s}$ can be $\texttt{[]}$.
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
