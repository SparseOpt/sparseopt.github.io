---
layout: archive
title: ""   
permalink: /SCO/
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

 
## Sparsity-constrained optimization
---

<p style="line-height: 1;"></p>
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n}} ~~  f(\mathbf{x}),~~~~ \mbox{s.t.}~~ \parallel\mathbf{x}\parallel_0\leq s  \tag{SCO}
\end{equation}

<div style="text-align:justify;">
where  $f:\mathbb{R}^{n}\rightarrow \mathbb{R}$ is a continuously or twice continuously differentiable function, $\mathrm{s\ll n}$ is a given integer, and $\|\mathbf{x}\|_0$ denotes the so-called L0 norm, which counts the number of nonzero entries in $\mathbf{x}$.
</div>
 

---
<div style="text-align:justify;"> 
The package can be downloaded here -<a style="font-size: 16px; font-weight: bold;color:#006DB0" href="\files\SCOpack.zip" target="_blank">SCOpack</a>,
which provides 3 solvers from the following papers:
</div>

> <b style="font-size:14px;color:#777777">NHTP</b> - <span style="font-size: 14px"> S Zhou, N Xiu, and H Qi, Global and quadratic convergence of Newton hard-thresholding pursuit, JMLR, 22:1-45, 2021. </span>
<br><b style="font-size:14px;color:#777777">GPNP</b> - <span style="font-size: 14px"> S Zhou, Gradient projection Newton pursuit for sparsity constrained optimization, ACHA, 61:75-100, 2022. </span>
<br><b style="font-size:14px;color:#777777">IIHT</b> - <span style="font-size: 14px"> L Pan, S Zhou, N Xiu, and H Qi, A convergent iterative hard thresholding for nonnegative sparsity optimization, PJO, 13:325-353, 2017. </span>

---
<div style="text-align:justify;">  
Note that $\texttt{NHTP}$ and $\texttt{GPNP}$ are two second-order methods, requiring the objective, gradient, and sub-Hessian of $f$, while $\texttt{IIHT}$ is a first-order method, requiring the objective and gradient of $f$. Below is an example of how to uniformly define functions for three solvers to solve a simple SCO problem, where input $\texttt{x}$ is the variable, string variable $\texttt{key}$ specifies the computation: $\texttt{key}$='$\texttt{fg}$' for the objective value and the gradient, and $\texttt{key}$='$\texttt{h}$' for the sub-Hessian decided by two indix sets $\texttt{T1}$ and $\texttt{T2}$. 

Specifically, when $\texttt{key}$ = '$\texttt{fg}$', if there is only one output, it returns the objective function value; if there are two outputs, the first is the objective function value and the second is the gradient. When $\texttt{key}$ = '$\texttt{h}$', if there is only one output, it returns the sub-Hessian consisting of the $\texttt{T1}$ rows and $\texttt{T1}$ columns of the Hessian; if there are two outputs, the first corresponds to the sub-Hessian formed by the $\texttt{T1}$ rows and $\texttt{T1}$ columns of the Hessian, and the second corresponds to the sub-Hessian formed by the $\texttt{T1}$ rows and $\texttt{T2}$ columns of the Hessian.
</div>
<p style="line-height: 1;"></p>

```ruby
function [out1,out2] = funcSimpleEx(x,key,T1,T2)
    % This code provides information for
    %     min   x'*[6 5;5 8]*x+[1 9]*x-sqrt(x'*x+1) 
    %     s.t. \|x\|_0<=s
    % where n=2 and s=1   
    a   = sqrt(sum(x.*x)+1);
    switch key
        case 'fg'    
            out1 = x'*[6 5;5 8]*x+[1 9]*x-a;       % objective
            if  nargout == 2 
                out2 = 2*[6 5;5 8]*x+[1; 9]-x./a;  % gradient
            end
        case 'h'
            H   = 2*[6 5;5 8]+(x*x'-a*eye(2))/a^3; % sub-Hessian formed by rows indexed by T1 and columns indexed by T1
            out1 = H(T1,T1);
            if  nargout == 2 
                out2 = H(T1,T2);                   % sub-Hessian formed by rows indexed by T1 and columns indexed by T2
            end
    end
end
```

<div style="text-align:justify;">
After defining the functions for the simple SCO problem, one can call $\texttt{SCOpack}$ to solve it. Users need to specify ($\texttt{func}$, $\texttt{n}$, $\texttt{s}$), choose a solver name from  {'$\texttt{NHTP}$', '$\texttt{GPNP}$', '$\texttt{IIHT}$'}, set some parameters in $\texttt{pars}$ if necessary, and then run the solver. The following codes demonstrate $\texttt{SCOpack}$ to solve this simple SCO problem.
</div>
<p style="line-height: 1;"></p>

```ruby
% demon a simple sparsity constrained problem
clc; close all; clear all;  addpath(genpath(pwd));

n        = 2;
s        = 1; 
func     = @funcSimpleEx;
solver   = {'NHTP','GPNP','IIHT'};
pars.eta = 0.1; % useful for 'NHTP'
out      = SCOpack(func,n,s,solver{2},pars); 

fprintf(' Objective:      %.4f\n', out.obj); 
fprintf(' CPU time:      %.3fsec\n', out.time);
fprintf(' Iterations:        %4d\n', out.iter);
```

<div style="text-align:justify;">
For other problems, users can similarly define the functions by modifying $\texttt{out1}$ and $\texttt{out2}$ while preserving the overall structure of $\texttt{switch}$. As an illustration, the following codes define the functions of a sparse linear regression problem.
</div>
<p style="line-height: 1;"></p>

```ruby
function [out1,out2] = funcLinReg(x,key,T1,T2,A,b)
    % This code provides information for
    %     min   0.5*||Ax-b||^2 
    %     s.t. \|x\|_0<=s
    % where A in R^{m x n} and b in R^{m x 1}    
    switch key
        case 'fg'
            Tx   = find(x~=0);
            Axb  = A(:,Tx)*x(Tx)-b;
            out1 = (Axb'*Axb)/2;      % objective 
            if  nargout == 2 
                out2 = (Axb'*A)';     % gradient 
            end
        case 'h'        
            AT   = A(:,T1); 
            out1 = AT'*AT;            %sub-Hessian formed by rows indexed by T1 and columns indexed by T1   
            if  nargout == 2
                out2 = AT'*A(:,T2);   %sub-Hessian formed by rows indexed by T1 and columns indexed by T2
            end       
    end
end
```
<div style="text-align:justify;">
After defining the functions of the sparse linear regression problem, we call $\texttt{SCOpack}$ to solve the problem as follows.
</div>
<p style="line-height: 1;"></p>

```ruby
% demon sparse linear regression problems 
clc; close all; clear all; addpath(genpath(pwd));

n        = 20000;  
m        = ceil(0.25*n); 
s        = ceil(0.025*n);

Tx       = randperm(n,s);  
xopt     = zeros(n,1);  
xopt(Tx) = randn(s,1); 
A        = randn(m,n)/sqrt(m); 
b        = A*xopt;  

func     = @(x,key,T1,T2)funcLinReg(x,key,T1,T2,A,b);
pars.tol = 1e-6;
solver   = {'NHTP','GPNP','IIHT'};
out      = SCOpack(func,n,s,solver{2},pars);
PlotRecovery(xopt,out.sol,[900,500,500,250],1)
```

<div style="text-align:justify;">
The inputs and outputs of $\texttt{SCOpack}$ are detailed below, where inputs ($\texttt{func}$, $\texttt{n}$, $\texttt{s}$, $\texttt{solvername}$) are required. The parameters in $\texttt{pars}$ are optional, but setting certain ones may improve the solver's performance and the solution quality. For example, if solver $\texttt{NHTP}$ is chosen, then tuning a proper $\texttt{pars.eta}$ may significantly improve the solver performance in terms of convergence speed and accuracy. Moreover, solver $\texttt{IIHT}$ enables addressing the SCO problems with non-negative constraints. To do so, just set $\texttt{pars.neg}$=1. 
</div>

<p style="line-height: 1;"></p>

```ruby
function out = SCOpack(func,n,s,solvername,pars)
% -------------------------------------------------------------------------
% This code aims at solving the sparsity constrained optimization (SCO),
%
%         min_{x\in R^n} f(x),  s.t. ||x||_0<=s
%
% or non-negative and sparsity constrained optimization (NSCO):
%
%         min_{x\in R^n} f(x),  s.t. ||x||_0<=s, x>=0 
%
% where f: R^n->R and s<<n is an integer.
% -------------------------------------------------------------------------
% Inputs:
%   func:   A function handle defines                            (REQUIRED)
%                    (objective, gradient, sub-Hessian)
%   n:      Dimension of the solution x                          (REQUIRED)
%   s:      Sparsity level of x, an integer between 1 and n-1    (REQUIRED)
%   solver: A text string, can be one of {'NHTP','GPNP','IIHT'}  (REQUIRED)
%   pars  : ---------------For all solvers --------------------------------
%           pars.x0    --  Starting point of x         (default zeros(n,1))
%           pars.disp  --  =1 show results for each step        (default 1)
%                          =0 not show results for each step
%           pars.maxit --  Maximum number of iterations      (default  2e3) 
%           pars.tol   --  Tolerance of halting conditions   (default 1e-6)
%           pars.uppf  --  An upper bound of final objective (default -Inf)
%                          Useful for noisy case
%           ---------------Particular for NHTP ----------------------------
%           pars.eta   --  A positive scalar                    (default 1)  
%                          Tuning it may improve solution quality 
%           ---------------Particular for IIHT ----------------------------
%           pars.neg   --  =0 for model (SCO)                   (default 1)
%                          =1 for model (NSCO)
% -------------------------------------------------------------------------
% Outputs:
%     out.sol :   The sparse solution x
%     out.obj :   Objective function value at out.sol 
%     out.iter:   Number of iterations
%     out.time:   CPU time
% -------------------------------------------------------------------------
% Send your comments and suggestions to <<< slzhou2021@163.com >>>   
% WARNING: Accuracy may not be guaranteed!!!!!  
% -------------------------------------------------------------------------
```
