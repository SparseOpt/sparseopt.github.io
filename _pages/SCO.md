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

 

##  <span style="color:#8C8C8C"> Sparsity-constrained optimization</span> 
---

<p style="line-height: 1;"></p>
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n}} ~~  f(\mathbf{x}),~~~~ \mbox{s.t.}~~ \parallel\mathbf{x}\parallel_0\leq s  \tag{SCO}
\end{equation}

<div style="text-align:justify;">
where  $f:\mathbb{R}^{n}\rightarrow \mathbb{R}$ is a continuously or twice continuously differentiable function, $\mathrm{s\ll n}$ is a given integer, and $\|\mathbf{x}\|_0$ denotes the so-called L0 norm, which counts the number of nonzero entries in $\mathbf{x}$.
</div>
 
<!-- ## <span style="color:#8C8C8C"> The solver and its demonstration </span> -->

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
Note that <b style="font-size:16px;color:#777777">NHTP</b> and <b style="font-size:16px;color:#777777">GPNP</b> are two second-order methods, requiring both the gradient and Hessian of $f$. Below is an example of how to define these for three solvers in the context of (SCO) problems uniformly, where input $\texttt{x}$ is the variable, string variable $\texttt{key}$ specifies the computation: $\texttt{key}$='$\texttt{fg}$' for the objective value and the gradient, and $\texttt{key}$='$\texttt{h}$' for the Hessian, and $\texttt{T1}$ and $\texttt{T2}$ are two indices and only valid when $\texttt{key}$='$\texttt{h}$'.
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
            H   = 2*[6 5;5 8]+(x*x'-a*eye(2))/a^3; % sub-Hessian on (T1 T1) 
            out1 = H(T1,T1);
            if  nargout == 2 
                out2 = H(T1,T2);                   % sub-Hessian on (T1 T2) 
            end
    end
end
```

<div style="text-align:justify;">
  Below is an example showing how <b style="font-size:16px;color:#777777">NM01</b> can be applied to solve the 1BCS problem using model (<a style="font-size: 16px;color:#006DB0" href="https://sparseopt.github.io/1BCS/" target="_blank">SFRO</a>). Users only need to specify ($\texttt{func}$, $\texttt{B}$, $\texttt{b}$, $\texttt{lam}$, $\texttt{pars}$) and then run the solver.
</div>

<p style="line-height: 1;"></p>

```ruby
% Solving 1-bit compressive sensing using randomly generated datasets
clc; close all; clear all; addpath(genpath(pwd));

n            = 1000; 
m            = ceil(0.5*n);
s            = ceil(0.01*m);                     % sparsity level
r            = 0.01;                             % flipping ratio
nf           = 0.1;                              % noisy ratio
[A,c,co,xo]  = random1bcs('Ind',m,n,s,nf,r,0.5); % data generation

func         = @(x,key)func1BCS(x,key,1e-5,0.5,A,c);
B            = (-c).*A;
b            = (n*4e-5)*ones(m,1);
lam          = 10;
pars.tau     = 1;
pars.sp      = s;  
pars.strict  =(n<=2000); 
out          = NM01(func, B, b, lam, pars); 
x            = refine(out.sol,s,A,c);

RecoverShow(xo,x,[950,500,450,220],1)
fprintf(' Computational time:    %.3fsec\n',out.time);
fprintf(' Signal-to-noise ratio: %.2f\n',-20*log10(norm(x-xo)));
fprintf(' Hamming distance:      %.3f\n',nnz(sign(A*x)-c)/m)
fprintf(' Hamming error:         %.3f\n',nnz(sign(A*x)-co)/m)
```

<div style="text-align:justify;">
The inputs and outputs of <b style="font-size:16px;color:#777777">NM01</b> are detailed below, where inputs ($\texttt{func}$, $\texttt{B}$, $\texttt{b}$, $\texttt{lam}$) are required. The parameters in $\texttt{pars}$ are optional, but setting certain ones may improve the solver's performance and the quality of the solution.
</div>

<p style="line-height: 1;"></p>

```ruby
function out = NM01(func,B,b,lam,pars)
% -------------------------------------------------------------------------
% This code aims at solving the support vector machine with form
%
%       min  f(x) + lam * ||(Bx+b)_+||_0
%
% where f is twice continuously differentiable
% lam > 0, B\in\R^{m x n}, b\in\R^{m x 1}
% (z)_+ = (max{0,z_1},...,max{0,z_m})^T
% ||(z)_+ ||_0 counts the number of positive entries of z
% -------------------------------------------------------------------------
% Inputs:
%   func: A function handle defines (objective,gradient,Hessain) (REQUIRED)
%   B   : A matrix \R^{m x n}                                    (REQUIRED)      
%   b   : A vector \R^{m x 1}                                    (REQUIRED)
%   lam : The penalty parameter                                  (REQUIRED)
%   pars: Parameters are all OPTIONAL
%         pars.x0     -- The initial point             (default zeros(n,1))
%         pars.tau    -- A useful paramter                   (default 1.00)
%         pars.mu0    -- A smoothing parameter               (default 0.01)
%         pars.maxit  -- Maximum number of iterations        (default 1000)  
%         pars.tol    -- Tolerance of halting conditions   (1e-7*sqrt(n*m)) 
%         pars.strict -- = 0, loosely meets halting conditions  (default 0)
%                        = 1, strictly meets halting conditions  
%                        pars.strict=1 is useful for low dimensions                           
% -------------------------------------------------------------------------
% Outputs:
%   out.sol:  The solution 
%   out.obj:  The objective function value
%   out.time: CPU time
%   out.iter: Number of iterations
% -------------------------------------------------------------------------
% Send your comments and suggestions to <<< slzhou2021@163.com >>>                                  
% WARNING: Accuracy may not be guaranteed!!!!!  
% -------------------------------------------------------------------------
```
