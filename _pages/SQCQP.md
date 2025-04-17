---
layout: archive
title: ""   
permalink: /SQCQP/
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

 

##  <span style="color:#8C8C8C"> Sparse quadratically constrained quadratic programming solver </span> 
---
Sparse quadratically constrained quadratic programming (SQCQP) takes the form of

\begin{equation}
\begin{aligned}
\min_{\mathbf{x}\in\mathbb{R}^{n}} &  \frac{1}{2}\mathbf{x}^{\top}\mathbf{Q}_0\mathbf{x}+\mathbf{q}_0^{\top}\mathbf{x}\\
\mbox{s.t.} & \frac{1}{2}\mathbf{x}^{\top}\mathbf{Q}_i\mathbf{x}+\mathbf{q}^{\top}_i\mathbf{x}+c_i\leq0,~i=1,2,\ldots,k\\
&\mathbf{A}\mathbf{x}\leq \mathbf{b}\\
&\mathbf{Aeq}\mathbf{x} = \mathbf{beq}\\
& lb\leq x_i \leq ub,~i=1,2,\ldots,n\\
&\parallel\mathbf{x}\parallel_0\leq s
\begin{aligned} \tag{SCO}
\end{equation}
where parameters are defined as follows
- $\mathbf{Q}_i\in\mathbb{R}^{n\times n}, \mathbf{q}_i\in\mathbb{R}^{n}, c_i\in\mathbb{R},~~i=0,1,\ldots,k$
- $\mathbf{A}\in\mathbb{R}^{m_1\times n}$, $\mathbf{b}\in\mathbb{R}^{m_1}$
- $\mathbf{Aeq}\in\mathbb{R}^{m_2\times n}$, $\mathbf{beq}\in\mathbb{R}^{m_2}$
- $lb$ and $ub$ are two scalars satisfying $0\in[lb, ub]$
- $\|\mathbf{x}\|_0$ denotes the so-called $\ell_0$-norm, which counts the number of nonzero entries in $\mathbf{x}$
- $s\ll n$ is a given integer
         
## <span style="color:#8C8C8C"> The solver and its demonstration </span> 
---
<div style="text-align:justify;">
<a style="font-size: 16px; font-weight: bold;color:#006DB0" href="https://github.com/ShenglongZhou/CSpack" target="_blank">SNSQP</a> was developed based on the algorithm developed in the following three papers:
</div>

<p style="line-height: 1;"></p>

- <a style="font-size: 14px;color:#000000" href="https://jmlr.org/papers/v22/19-026.html" target="_blank"> S Zhou, N Xiu and H  Qi, Global and quadratic convergence of Newton hard-thresholding pursuit, *J Mach Learn Res*, 22:1−45, 2021.</a>
- <a style="font-size: 14px;color:#000000" href="https://www.sciencedirect.com/science/article/pii/S1063520322000458" target="_blank"> S Zhou, Gradient projection newton pursuit for sparsity constrained optimization, *Appl Comput Harmon Anal*, 61:75-100, 2022.</a> 
- <a style="font-size: 14px;color:#000000" href="http://www.yokohamapublishers.jp/online2/oppjo/vol13/p325.html" target="_blank"> L Pan, S Zhou, N Xiu, and H Qi, A convergent iterative hard thresholding for nonnegative sparsity optimization, *Pac J Optim*, 13:325-353, 2017.</a>
 

<b style="font-size:14px;color:#777777">NHTP</b> and <b style="font-size:14px;color:#777777">GPNP</b> are second-order methods, which require both the gradient and Hessian of $f$. In contrast, <b style="font-size:14px;color:#777777">IIHT</b> is a first-order method that only requires the gradient of $f$. Below is a demonstration of how to define the gradient and Hessian for <b style="font-size:14px;color:#777777">NHTP</b>.

<p style="line-height: 1;"></p>

```ruby
function [out1,out2] = funCS(x,T1,T2,data)

    if  isempty(T1) && isempty(T2) 
        Tx   = find(x); 
        Axb  = data.A(:,Tx)*x(Tx)-data.b;
        out1 = norm(Axb,'fro')^2/2;               %objective 
        if  nargout == 2
            out2    = (Axb'*data.A)';             %gradient
        end
    else        
        AT = data.A(:,T1); 
        if  length(T1)<2000
            out1 = AT'*AT;                        %subHessian containing T1 rows and T1 columns
        else
            out1 = @(v)( (AT*v)'*AT )';      
        end       
        if  nargout == 2
            out2 = @(v)( (data.A(:,T2)*v)'*AT )'; %subHessian containing T1 rows and T2 columns
        end       
    end     
end
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
 
