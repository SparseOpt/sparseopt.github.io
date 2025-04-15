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

 

##  <span style="color:##555555"><b> Sparsity-constrained optimization solvers</b></span> 
---
Sparsity-constrained optimization (SCO) takes the form of

\begin{eqnarray}\min_{\mathbf{x}} &&   f(\mathbf{x}) \nonumber\\\\\\
\mbox{s.t.} && \parallel\mathbf{x}\parallel_0\leq s,~\mathbf{x}\in\mathbb{R}^{n}, \nonumber
\end{eqnarray}
<div style="text-align:justify;">
where  $f:\mathbb{R}^{n}\rightarrow \mathbb{R}$ is a continuously or twice continuously differentiable function, $s\ll n$ is a given integer, and $\|\mathbf{x}\|_0$ denotes the so-called $\ell_0$-norm, which counts the number of nonzero entries in $\mathbf{x}$.
</div>
 
### <span style="color:##555555"><b style="font-size:20px">The solver and its demonstration</b></span> 
---

<div style="text-align:justify;">
<a style="font-size: 16px; font-weight: bold;color:#006DB0" href="https://github.com/ShenglongZhou/CSpack" target="_blank">SCOsolvers</a> provides three solvers:  <b style="font-size:16px;color:#555555">NHTP</b>, <b style="font-size:16px;color:#555555">GPNP</b>, and  <b style="font-size:16px;color:#555555">IIHT</b> based on the algorithms developed in the following three papers:
</div>

<p style="line-height: 1;"></p>

```
S Zhou, N Xiu and H  Qi, Global and quadratic convergence of Newton hard-thresholding pursuit, J Mach Learn Res, 22(12):1−45, 2021. 

S Zhou, Gradient projection newton pursuit for sparsity constrained optimization, Appl Comput Harmon Anal,  61:75-100, 2022. 

L Pan, S Zhou, N Xiu, and H Qi, A convergent iterative hard thresholding for nonnegative sparsity optimization, Pac J Optim, 13(2):325-353, 2017.
```

<b style="font-size:14px;color:#555555">NHTP</b> and <b style="font-size:14px;color:#555555">GPNP</b> are second-order methods, which require both the gradient and Hessian of $f$. In contrast, <b style="font-size:14px;color:#555555">IIHT</b> is a first-order method that only requires the gradient of $f$. Below is a demonstration of how to define the gradient and Hessian for <b style="font-size:14px;color:#555555">NHTP</b> and <b style="font-size:14px;color:#555555">GPNP</b>.

<p style="line-height: 1;"></p>

```ruby
clc; clear; close all; 

ExName     = 'DempeDutta2012Ex24'; 
func       = str2func(ExName);
dim        = [1 1 0 1];     % required
pars.xy    = [1;1];         % optional

pars.lam   = 1;             % optional
pars.keep  = 0;             % optional 
pars.check = 1;             % optional

SolNo      = 1;     % choose a solver
Solvers    = {'SNLLVF','SNQVI','SNKKT'}; 
solver     = str2func(Solvers{SolNo});  
Out1       = solver(func, dim,  pars);
```

<div style="text-align:justify;">
Each solver has two required inputs: 'func' defining the example and 'dim' recording dimensions of the example and an optional input 'pars'. Please see the <a style="font-size: 16px; font-weight: bold; color:#007D98" href="\files\menu-of-BiOpt.pdf" target="_blank">menu-of-BiOpt</a> for more details of 'pars'. The chosen solver is <span style="color:#007D98"><b style="font-size:16px">SNLLVF</b></span> and solves one example 'DempeDutta2012Ex24' defined by the following Matlab m-file. This example is from <a style="font-size: 16px; font-weight: bold; color:#007D98"  href="https://biopt.github.io/bolib/" target="_blank">BOLIBver2</a>, where more examples are provided. The <a style="font-size: 16px; font-weight: bold; color:#007D98" href="\files\menu-of-BiOpt.pdf" target="_blank">menu-of-BiOpt</a> also presents other ways to define examples.
</div>

<p style="line-height: 1;"></p>

```ruby
function w=DempeDutta2012Ex24(x,y,keyf,keyxy)
% This file provides all functions defining DempeDutta2012Ex24 problem and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 0 1]
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = (x-1)^2+y^2;
    case 'G'; w = []; 
    case 'f'; w = x^2*y;      
    case 'g'; w = y^2; 
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = 2*(x-1);         
        case 'y' ; w = 2*y;        
        case 'xx'; w = 2;
        case 'xy'; w = 0;
        case 'yy'; w = 2;
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [];    
        case 'y' ; w = [];      
        case 'xx'; w = [];
        case 'xy'; w = [];
        case 'yy'; w = [];
        end           
    case 'f'   
        switch keyxy
        case 'x' ; w = 2*x*y;    
        case 'y' ; w = x^2;          
        case 'xx'; w = 2*y;
        case 'xy'; w = 2*x;
        case 'yy'; w = 0;
        end           
    case 'g'   
        switch keyxy
        case 'x' ; w =   0;  
        case 'y' ; w =   2*y;         
        case 'xx'; w =   0;  
        case 'xy'; w =   0;  
        case 'yy'; w =   2; 
        end        
   end   
end
end
```
 
