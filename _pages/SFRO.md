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
\min_{\mathbf{x}\in\mathbb{R}^{n}} ~~  f(\mathbf{x}) + \lambda \parallel(\mathbf{A}\mathbf{x}+\mathbf{b})_+\parallel_0  \tag{SFRO}
\end{equation}

<div style="text-align:justify;">
where  $f:\mathbb{R}^{n}\rightarrow \mathbb{R}$ is a twice continuously differentiable function,  $\mathbf{A}\in\mathbb{R}^{m\times n}$ is a matrix, $\mathbf{b}\in\mathbb{R}^{m}$ is a vector, $\lambda$ is a given scalar,  $\mathbf{z}_+=(\max\{0,z_1\},\ldots,\max\{0,z_m\})^\top$, and  $\|\mathbf{z}\|_0$ is the so-called $\ell_0$ norm that counts the number of nonzero entries in $\mathbf{z}$. The regularization term is related to the step function (or 0/1 loss function) defined by $\ell_{0/1}(t)=1$ if $t>0$ and $\ell_{0/1}(t)=0$ otherwise. Therefore,
  \begin{equation}\|\mathbf{z}_+\|_0=\sum_{i=1}^m \ell_{0/1}\left(z_{i}\right)\nonumber\end{equation}
</div>
 
<!-- ## <span style="color:#8C8C8C"> The solver and its demonstration </span> -->

---
<div style="text-align:justify;"> 
The solver can be downloaded here - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href=" " target="_blank">NM01</a>,
which was developed from the following paper:
</div>

> <span style="font-size: 14px"> S Zhou, L Pan, N Xiu,  and H Qi, Quadratic convergence of smoothing Newton's method for 0/1 loss optimization, SIOPT, 31:3184–3211, 2021. </span>


---
<div style="text-align:justify;">  
Note that <b style="font-size:14px;color:#777777">NM01</b> is a second-order method, which requires the gradient and Hessian of $f$. Below is a demonstration of how to define the gradient and Hessian for the solver.
</div>

```ruby
function out = funcSVM(x,key,w,A0,c)
    n = length(x);
    switch key   
        case 'f';  out = norm(x,'fro')^2 - (1-w)*x(n)^2;
        case 'g';  out = x;  out(n) = w*x(n);
        case 'h';  out = speye(n); out(n,n) = w; 
        case 'a';  acc = @(var)nnz( sign( (A0*var(1:n-1)+var(n)) )-c);
                   out = 1-acc(x)/length(c);
        otherwise; out = []; % 'Otherwise' is REQIURED
    end    
end
```

---
<div style="text-align:justify;">
Below is a demonstration of how <b style="font-size:14px;color:#777777">NM01</b> can be used to solve the problem. You simply need to input the data $(\texttt{A},\texttt{b},\texttt{s},\texttt{k})$ and then choose one solver from $\texttt{\{`GPSP',`NM01'\}}$. 
</div>

<p style="line-height: 1;"></p>

```ruby
% Solving support vector machine using four synthetic samples
clc; close all; clear all;  addpath(genpath(pwd));

a          = 10;
A0          = [0 0; 0 1; 1 0; 1 a]; 
c           = [-1 -1  1  1]';
[m,n]       = size(A0);  

func        = @(x,key)funcSVM(x,key,1e-4,A0,c);
A           = (-c).*[A0 ones(m,1)];
b           = ones(m,1);
lam         = 10;
pars.tau    = 1;
pars.strict = 1;
out         = NM01(func, A, b, lam, pars); 
x           = out.sol;        

figure('Renderer', 'painters', 'Position', [1000, 300,350 330])
axes('Position', [0.08 0.08 0.88 0.88] );
scatter([1;1],[0 a],80,'+','m'), hold on
scatter([0;0],[0,1],80,'x','b'), hold on
line([-x(3)/x(1) -x(3)/x(1)],[-1 1.1*a],'Color', 'r')
axis([-.1 1.1 -1 1.1*a]),box on
ld = strcat('NM01:',num2str(func(x,'a')*100,'%.0f%%'));
legend('Positive','Negative',ld,'location','NorthWest')
```

<div style="text-align:justify;">
The inputs and outputs of NM01 are detailed below, where inputs $(\texttt{A},\texttt{b},\texttt{s},\texttt{k},\texttt{solver})$ are required. If choose $\texttt{solver=`NM01'}$, then one can set $\texttt{s=[]}$ and $ \texttt{k=[]}$ if they are unknown. The parameters in $\texttt{pars}$ are optional, but setting certain ones may improve the solver's performance and the quality of the solution.
</div>

<p style="line-height: 1;"></p>

```ruby
function out = NM01(func,A,b,lam,pars)
% -------------------------------------------------------------------------
% This code aims at solving the support vector machine with form
%
%      min  f（x） + lam * ||(Ax+b)_+||_0
%
% where f is twice continuously differentiable
% lam > 0, A\in\R^{m x n}, b\in\R^{m x 1}
% (z)_+ = (max{0,z_1},...,max{0,z_m})^T
% ||(z)_+ ||_0 counts the number of positive entries of z
% -------------------------------------------------------------------------
% Inputs:
%   func: A function handle defines (objective,gradient,Hessain) (REQUIRED)
%   A   : A matrix \R^{m x n}                                    (REQUIRED)      
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
