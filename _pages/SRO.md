---
layout: archive
title: ""   
permalink: /SRO/
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


## Sparsity-regularized optimization
---
<p style="line-height: 2;"></p>

\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n}} ~~  f(\mathbf{x}) + \lambda \parallel\mathbf{x}\parallel_0 \tag{SRO}
\end{equation}

<div style="text-align:justify;"> 
where $f:\mathbb{R}^{n}\rightarrow \mathbb{R}$ is a continuously or twice continuously differentiable function, $\lambda>0$ is the penalty parameter, and $\parallel\mathbf{x}\parallel_0$ denotes the L0 norm, which counts the number of nonzero entries in $\mathbf{x}$.
</div>
 
---
<div style="text-align:justify;">
Package - <a style="font-size: 16px; font-weight: bold; color:#006DB0" href="\files\SROpack-Matlab.zip" target="_blank">SROpack-Matlab</a> and <a style="font-size: 16px; font-weight: bold; color:#006DB0" href="\files\SROpack-Python.zip" target="_blank">SROpack-Python</a> (click to download) provides 1 solver from the following paper: 
</div>

> <b style="font-size:14px;color:#777777">NL0R</b> - <span style="font-size: 14px"> S Zhou, L Pan, and N Xiu, Newton method for l0 regularized optimization, Numer Algorithms, 88:1541–1570, 2021. </span>

---
<div style="text-align:justify;">  
Note that $\texttt{NL0R}$ is a second-order method, requiring the objective, gradient, and sub-Hessian of $f$. Based on Matlab syntax  (similar to Python syntax), below is an example of how to define these functions for the solver to solve a simple problem, where input $\texttt{x}$ is the variable, $\texttt{key}$ is a string variable, and $\texttt{T1}$ and $\texttt{T2}$ are two index sets.  Here, $\texttt{key}$ is used to specify the computation: when $\texttt{key}$ = '$\texttt{f}$', output the objective function value; when $\texttt{key}$ = '$\texttt{g}$', output the gradient; when $\texttt{key}$ = '$\texttt{h}$', output the sub-Hessian containing the $\texttt{T1}$ rows and $\texttt{T2}$ columns of the Hessian. 
</div>
<p style="line-height: 1;"></p>

```ruby
function  out = funcSimpleEx(x,key,T1,T2)
    % This code provides information for
    %     min   x'*[6 5;5 8]*x+[1 9]*x-sqrt(x'*x+1) 
    a   = sqrt(sum(x.*x)+1);
    switch key
        case 'f'    
            out = x'*[6 5;5 8]*x+[1 9]*x-a;         % objective
        case 'g'    
            out = 2*[6 5;5 8]*x+[1; 9]-x./a;        % gradient
        case 'h'
            H   = 2*[6 5;5 8]+(x*x'-a*eye(2))/a^3;  % sub-Hessian indexed by T1 and T2 
            out = H(T1,T2);
    end
end
```

<div style="text-align:justify;">
After defining the functions for the simple problem, one can call $\texttt{NL0R}$ to solve it. Users need to specify ($\texttt{func}$, $\texttt{n}$, $\texttt{lambda}$), set some parameters in $\texttt{pars}$ if necessary, and then run the solver. The following Matlab code demonstrates $\texttt{NL0R}$ to solve this simple problem.
</div>
<p style="line-height: 1;"></p>

```ruby
% demon a simple SRO problem
clc; close all; clear all;  addpath(genpath(pwd));

n        = 2;
lambda   = 0.5;
pars.eta = 0.1;
out      = NL0R(@funcSimpleEx,n,lambda,pars); 
fprintf(' Objective:      %.4f\n', out.obj); 
fprintf(' CPU time:      %.3fsec\n', out.time);
fprintf(' Iterations:        %4d\n', out.iter);
```

<div style="text-align:justify;">
For other problems, users can similarly define the functions by modifying $\texttt{out1}$ and $\texttt{out2}$ while preserving the overall structure of $\texttt{switch}$. For example, the following Matlab code defines functions for a sparse linear regression problem. Similarly, in the inputs of function $\texttt{funcLinReg}$, ($\texttt{x}$, $\texttt{key}$, $\texttt{T1}$, $\texttt{T2}$) are variables, while ($\texttt{A}$, $\texttt{b}$) are data. When calling function $\texttt{funcLinReg}$, users need to define these data.
</div>
<p style="line-height: 1;"></p>

```ruby
function out = funcLinReg(x,key,T1,T2,A,b)
    % This code provides information for
    %     min   0.5*||Ax-b||^2 
    % where A in R^{m x n} and b in R^{m x 1}    
    switch key
        case 'f'
            Tx   = find(x~=0);
            Axb  = A(:,Tx)*x(Tx)-b;
            out  = (Axb'*Axb)/2;             % objective  
        case 'g'
            Tx   = find(x~=0); 
            out  = ((A(:,Tx)*x(Tx)-b)'*A)';  % gradient   
        case 'h'        
            out  = A(:,T1)'*A(:,T2);         % sub-Hessian indexed by T1 and T2
end    
```

<div style="text-align:justify;">
After defining the functions of the sparse linear regression problem, we call $\texttt{NL0R}$ to solve the problem as follows.
</div>
<p style="line-height: 1;"></p>

```ruby
% demon sparse linear regression problems 
clc; close all; clear all; addpath(genpath(pwd));

n        = 2000;  
m        = ceil(0.25*n); 
s        = ceil(0.05*n);

Tx       = randperm(n,s);  
xopt     = zeros(n,1);  
xopt(Tx) = (0.25+rand(s,1)).*sign(randn(s,1)); 
A        = randn(m,n)/sqrt(m); 
b        = A*xopt;  
func     = @(x,key,T1,T2)funcLinReg(x,key,T1,T2,A,b);

lambda   = 0.01;
pars.eta = 1.0;
out      = NL0R(func,n,lambda,pars); 
```

<div style="text-align:justify;">
The inputs and outputs of the Matlab version of $\texttt{NL0R}$ are detailed below, with analogous specifications for the Python version. Inputs ($\texttt{func}$, $\texttt{n}$, $\texttt{lambda}$) are required. The parameters in $\texttt{pars}$ are optional, but setting certain ones may improve the solver's performance and the solution quality. For example, tuning a proper $\texttt{pars.eta}$ may significantly improve the solver performance in terms of convergence speed and accuracy.
</div>

<p style="line-height: 1;"></p>

```ruby
function out = NL0R(func,n,lambda,pars)
% -------------------------------------------------------------------------
% This code aims at solving the L0 norm regularized optimization 
%
%         min_{x\in R^n} f(x) + lambda*||x||_0^0
%
% where f: R^n->R, lambda>0
% ||x||_0^0 counts the number of non-zero entries
% -------------------------------------------------------------------------
% Inputs:
%   func:   A function handle defines                            (REQUIRED)
%                    (objective, gradient, sub-Hessian)
%   n:      Dimension of the solution x                          (REQUIRED) 
%   lambda: The penalty parameter                                (REQUIRED)  
%   pars  : All parameters are OPTIONAL
%           pars.x0     -- Starting point of x         (default zeros(n,1))
%           pars.tol    -- Tolerance of halting conditions   (default 1e-6)
%           pars.maxit  -- Maximum number of iterations      (default  2e3) 
%           pars.uppf   -- An upper bound of final objective (default -Inf)
%                          Useful for noisy case 
%           pars.eta    -- A positive scalar                    (default 1)  
%                          Tuning it may improve solution quality
%           pars.update -- =1 update penalty parameter lambda   (default 1)
%                          =0 fix penalty parameter lambda
%           pars.disp   -- =1 show results for each step        (default 1)
%                          =0 not show results for each step
% -------------------------------------------------------------------------
% Outputs:
%   out.sol :  The sparse solution x
%   out.obj :  Objective function value at out.sol 
%   out.iter:  Number of iterations
%   out.time:  CPU time
% -------------------------------------------------------------------------
% Send your comments and suggestions to <<< slzhou2021@163.com >>>   
% WARNING: Accuracy may not be guaranteed!!!!!  
% -------------------------------------------------------------------------
```
