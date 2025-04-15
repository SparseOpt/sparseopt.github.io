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

 

### <span style="color:#8C8C8C"><b style="font-size:20px">Sparsity-regularized optimization solvers</b></span> 
---
Sparsity-regularized optimization (SRO) takes the form of

\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n}} ~~  f(\mathbf{x}) + \lambda \parallel\mathbf{x}\parallel_q^q \tag{SRO}
\end{equation}

<div style="text-align:justify;">
where  $f:\mathbb{R}^{n}\rightarrow \mathbb{R}$ is a continuously or twice continuously differentiable function, $\lambda>0$ is a given scalar, and $\|\mathbf{x}\|_q^q=\sum_i |x_i|^q$ with $q\in[0,1)$ denotes the $\ell_q$-norm. In particular, when $q=0$,  $\|\mathbf{x}\|_0=\|\mathbf{x}\|_0^0$ is the so-called $\ell_0$-norm that counts the number of nonzero entries in $\mathbf{x}$.
</div>
 
### <span style="color:#8C8C8C"><b style="font-size:20px">The solver and its demonstration</b></span> 
---

<div style="text-align:justify;">
<a style="font-size: 16px; font-weight: bold; color:#006DB0" href="https://github.com/ShenglongZhou/CSpack" target="_blank">SROsolvers</a> provides two solvers:  <b style="font-size:14px;color:#777777">NL0R</b> and <b style="font-size:14px;color:#777777">PSNP</b> based on the algorithms developed in the following two papers:
</div>

- <a style="font-size:14px; color:#000000" href="https://link.springer.com/article/10.1007/s11075-021-01085-x" target="_blank"> S Zhou, L Pan, and N Xiu, Newton method for L0-regularized optimization, *Numerical Algorithm*, 88:1541–1570, 2021 .</a>
- <a style="font-size:14px; color:#000000" href="https://arxiv.org/abs/2306.14394" target="_blank"> S Zhou, X Xiu, Y Wang, and D Peng, Revisiting Lq (0 <= q < 1) norm regularized optimization, *arXiv:2306.14394*, 2023.</a>



<p style="line-height: 1;"></p>

Both solvers are second-order methods, which require both the gradient and Hessian of $f$. Below is a demonstration of how to define the gradient and Hessian for <b style="font-size:14px;color:#777777">NL0R</b> and <b style="font-size:14px;color:#777777">PSNP</b>.

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
 
