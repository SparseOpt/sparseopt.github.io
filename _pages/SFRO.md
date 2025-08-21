---
layout: archive
title: ""   
permalink: /SFRO/
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

 

##  <span style="color:#8C8C8C"> Step function-regularized optimization</span> 
---

<p style="line-height: 1;"></p>
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n}} ~~  f(\mathbf{x}) + \lambda \parallel(\mathbf{B}\mathbf{x}+\mathbf{b})_+\parallel_0  \tag{SFRO}
\end{equation}

<div style="text-align:justify;">
where  $f:\mathbb{R}^{n}\rightarrow \mathbb{R}$ is a twice continuously differentiable function,  $\mathbf{B}\in\mathbb{R}^{m\times n}$ is a matrix, $\mathbf{b}\in\mathbb{R}^{m}$ is a vector, $\lambda$ is a penalty parameter,  $\mathbf{z}_+=(\max\{0,z_1\},\ldots,\max\{0,z_m\})^\top$, and  $\|\mathbf{z}\|_0$ is the L0 norm that counts the number of nonzero entries in $\mathbf{z}$. Therefore, $\|\mathbf{z}_+\|_0$ counts the number of positive entries in $\mathbf{z}$.  This term is related to the step (or 0/1 loss) function defined by $\mathrm{step}(t)=1$ if $t>0$ and $\mathrm{step}(t)=0$ otherwise. As a result,  $\|\mathbf{z}_+\|_0= \mathrm{step}(z_1)+\cdots+\mathrm{step}(z_m).$
</div>
 
<!-- ## <span style="color:#8C8C8C"> The solver and its demonstration </span> -->

---
<div style="text-align:justify;"> 
The solver can be downloaded here - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href="\files\SFROpack.zip" target="_blank">NM01</a>,
which was developed from the following paper:
</div>

> <span style="font-size: 14px"> S Zhou, L Pan, N Xiu,  and H Qi, Quadratic convergence of smoothing Newton's method for 0/1 loss optimization, SIOPT, 31:3184–3211, 2021. </span>

 
---
<div style="text-align:justify;">  
Note that <b style="font-size:16px;color:#777777">NM01</b> is a second-order method, requiring both the gradient and Hessian of $f$. Below is an example of how to define these for the solver in the context of the 1-bit compressive sensing (<a style="font-size: 16px; font-weight: bold; color:#006DB0" href="https://sparseopt.github.io/1BCS/" target="_blank">1BCS</a>) problem, where objective function  $f(\mathbf{x})$ is given in model (<a style="font-size: 16px;color:#006DB0" href="https://sparseopt.github.io/1BCS/" target="_blank">SFRO</a>). The MATLAB codes below define $(f(\mathbf{x}), \nabla f(\mathbf{x}), \nabla^2 f(\mathbf{x}))$, where $\texttt{x}$ and $\texttt{key}$ are two variables, and  $\texttt{eps}$, $\texttt{q}$, $\texttt{A}$, and $\texttt{c}$ are parameters and data, as shown in model (<a style="font-size: 16px;color:#006DB0" href="https://sparseopt.github.io/1BCS/" target="_blank">SFRO</a>).  String variable $\texttt{key}$ specifies the computation: $\texttt{key}$='$\texttt{f}$' for the objective value, 
$\texttt{key}$='$\texttt{g}$' for the gradient, and $\texttt{key}$='$\texttt{h}$' for the Hessian. When $\texttt{key}$='$\texttt{a}$', an additional user-defined function is evaluated. Here,  the accuracy is computed for the 1BCS problem. This allows users to monitor a customized metric during optimization.
</div>
<p style="line-height: 1;"></p>

```ruby
function out = func1BCS(x,key,eps,q,A,c) 
    switch key   
        case 'f';  out = sum((x.^2+eps).^(q/2));
        case 'g';  out = q*x.*(x.^2+eps).^(q/2-1); 
        case 'h';  x2  = x.*x; out = diag(((x2+eps).^(q/2-2)).*((q-1)*x2+eps)); 
        case 'a';  acc = @(var)nnz(sign(A*var)-c); out = 1-acc(x)/length(c);
        otherwise; out = []; % 'Otherwise' is REQIURED if no key='a'
    end    
end
```

<div style="text-align:justify;">  
If no additional function is required, users can simply define $(f(\mathbf{x})$, $\nabla f(\mathbf{x})$, $\nabla^2 f(\mathbf{x}))$ by omitting case $\texttt{key}$='$\texttt{a}$' as follows.
</div>
<p style="line-height: 1;"></p>

```ruby
function out = func1BCS(x,key,eps,q) 
    switch key   
        case 'f';  out = sum((x.^2+eps).^(q/2));
        case 'g';  out = q*x.*(x.^2+eps).^(q/2-1); 
        case 'h';  x2  = x.*x; out = diag(((x2+eps).^(q/2-2)).*((q-1)*x2+eps)); 
        otherwise; out = [];  
    end    
end
```

<div style="text-align:justify;">
  Below is an example showing how <b style="font-size:16px;color:#777777">NM01</b> can be applied to solve the 1BCS problem using model (<a style="font-size: 16px;color:#006DB0" href="https://sparseopt.github.io/1BCS/" target="_blank">SFRO</a>). Users only need to specify ($\texttt{func}$, $\texttt{B}$, $\texttt{b}$, $\texttt{lam}$, $\texttt{pars}$) and then run the solver.
</div>

<p style="line-height: 1;"></p>

```ruby
% Solving 1 bit compressive sensing using randomly generated data 
clc; close all; clear all; addpath(genpath(pwd));

n            = 1000; 
m            = ceil(0.5*n);
s            = ceil(0.01*n);                      % sparsity level
r            = 0.01;                              % flipping ratio
nf           = 0.05;                              % noisy ratio
[A,c,co,xo]  = random1bcs('Ind',m,n,s,nf,r,0.5);  % data generation

func         = @(x,key)func1BCS(x,key,1e-5,0.5,A,c);
B            = (-c).*A;
b            = (n*8e-5)*ones(m,1);
lam          = 1;
pars.tau     = 1;  
pars.strict  =(n<=2000); 
out          = NM01(func, B, b, lam, pars); 
x            = refine(out.sol,s,A,c);

PlotRecovery(xo,x,[950,500,500,250],1)
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
