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

\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n}} ~~  f(\mathbf{x}),~~~~ \mbox{s.t.}~~ \parallel\mathbf{x}\parallel_0\leq s  \tag{SCO}
\end{equation}

<div style="text-align:justify;">
where  $f:\mathbb{R}^{n}\rightarrow \mathbb{R}$ is a continuously or twice continuously differentiable function, $s\ll n$ is a given integer, and $\|\mathbf{x}\|_0$ denotes the so-called $\ell_0$-norm, which counts the number of nonzero entries in $\mathbf{x}$.
</div>
 
<!-- ## <span style="color:#8C8C8C"> The solver and its demonstration </span> -->

---
<div style="text-align:justify;"> 
The package can be download here - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href="https://github.com/ShenglongZhou/CSpack" target="_blank">SCOsolvers</a>,
which provides 3 solvers from the following papers:
</div>

> <b style="font-size:14px;color:#777777">NHTP</b> - <span style="font-size: 14px"> S Zhou, N Xiu, and H Qi, Global and quadratic convergence of Newton hard-thresholding pursuit, J Mach Learn Res, 22:1-45, 2021. </span>
<br><b style="font-size:14px;color:#777777">GPSP</b> - <span style="font-size: 14px"> S Zhou, Gradient projection Newton pursuit for sparsity constrained optimization, Appl Comput Harmon Ana, 61:75-100, 2022. </span>
<br><b style="font-size:14px;color:#777777">IIHT</b> - <span style="font-size: 14px"> L Pan, S Zhou, N Xiu, and H Qi, A convergent iterative hard thresholding for nonnegative sparsity optimization, Pac J Optim, 13:325-353, 2017. </span>

<!--
- <a style="font-size: 14px;color:#000000" href="https://jmlr.org/papers/v22/19-026.html" target="_blank"> S Zhou, N Xiu and H  Qi, Global and quadratic convergence of Newton hard-thresholding pursuit, *J Mach Learn Res*, 22:1−45, 2021.</a>
- <a style="font-size: 14px;color:#000000" href="https://www.sciencedirect.com/science/article/pii/S1063520322000458" target="_blank"> S Zhou, Gradient projection newton pursuit for sparsity constrained optimization, *Appl Comput Harmon Anal*, 61:75-100, 2022.</a> 
- <a style="font-size: 14px;color:#000000" href="http://www.yokohamapublishers.jp/online2/oppjo/vol13/p325.html" target="_blank"> L Pan, S Zhou, N Xiu, and H Qi, A convergent iterative hard thresholding for nonnegative sparsity optimization, *Pac J Optim*, 13:325-353, 2017.</a> -->
 
<div style="text-align:justify;">  
Note that <b style="font-size:14px;color:#777777">NHTP</b> and <b style="font-size:14px;color:#777777">GPNP</b> are second-order methods, which require the gradient and Hessian of $f$. <b style="font-size:14px;color:#777777">IIHT</b> is a first-order method that only requires the gradient. Below is a demonstration of how to define the gradient and Hessian for <b style="font-size:14px;color:#777777">NHTP</b>.
</div>

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
