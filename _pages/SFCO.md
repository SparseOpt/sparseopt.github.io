---
layout: archive
title: ""   
permalink: /SFCO/
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

## Step function-constrained optimization
---

<p style="line-height: 1;"></p>
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{K}} ~~  f(\mathbf{x}),~~~~ \mbox{s.t.}~~ \parallel\mathbf{G}(\mathbf{x})\parallel_0^+\leq s,~\mathbf{x}\in \Omega  \tag{SFCO}
\end{equation}

<div style="text-align:justify;">
where $\mathbf{G}(\mathbf{x})\in\mathbb{R}^{M \times N}$ with each the ($i,j$)th entry being $G_{ij}(\mathbf{x})$, functions $f:\mathbb{R}^{K}\rightarrow \mathbb{R}$ and $G_{ij}:\mathbb{R}^{K}\rightarrow \mathbb{R}$ are (preferably twice) continuously differentiable, $s\ll N$ is an integer, and $\Omega\subseteq\mathbb{R}^{K}$ is a closed and convex set. In particular, it can be one of the following sets:
</div>
  
- A non-negative cone: $\{\mathbf{x}: x_i \geq 0, i=1,\ldots,K\}$ 
- A box with a lower bound $l$ and upper bound $u$:  $\{\mathbf{x}: l\leq x_i \leq u, i=1,\ldots,K\}$ 
- A sphare with a radius $r>0$: $\{\mathbf{x}: \|\mathbf{x}\|^2\leq r\}$ 
- A halfspace: $\{\mathbf{x}: \mathbf{a}^T\mathbf{x}\leq b\}$, where $\mathbf{a}\in\mathbb{R}^{K}$ and $b\in\mathbb{R}$
- A hyperplane: $\{\mathbf{x}: \mathbf{A} \mathbf{x}=  \mathbf{b}\}$, where $\mathbf{A}\in\mathbb{R}^{S\times K}$ and $ \mathbf{b}\in\mathbb{R}^{S}$

<div style="text-align:justify;">
  For a matrix $\mathbf{Z}\in\mathbb{R}^{M \times N}$,  measure $\|\mathbf{Z}\|_0^+$ counts the number of its columns that have positive entries, namely, 
  \begin{equation}\|\mathbf{Z}\|_0^+= \mathrm{step}\Big(\max_{i=1,\ldots,M} Z_{i1}\Big)+\cdots+\mathrm{step}\Big(\max_{i=1,\ldots,M} Z_{iN}\Big)\nonumber\end{equation}
  Here, $\mathrm{step}(t)$ is the step function (or 0/1 loss function) defined by $\mathrm{step}(t)=1$ if $t>0$ and $\mathrm{step}(t)=0$ otherwise. When $M=1$, matrix $\mathbf{Z}$ reduces to a vector $\mathbf{z}\in\mathbb{R}^{N}$, let $\mathbf{z}_+=$ $(\max\{0,z_1\}$ $\ldots$ $\max\{0,z_N\})^\top$ and  $\parallel\mathbf{z}\parallel_0$ denote its L0 norm that counts the number of its nonzero entries. As a result,  
  \begin{equation*}\|\mathbf{z}\|_0^+= \mathrm{step}(z_1)+\cdots+\mathrm{step}(z_N)=\|\mathbf{z}_+\|_0\end{equation*}
</div>
 
<!-- ## <span style="color:#8C8C8C"> The solver and its demonstration </span> -->

---
<div style="text-align:justify;"> 
Package - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href=" " target="_blank">SFCOpack-Matlab</a> or <a style="font-size: 16px; font-weight: bold;color:#006DB0" href=" " target="_blank">SFCOpack-Python</a>  provides 1 solver from the following paper:
</div>

> <b style="font-size:14px;color:#777777">SNSCO</b> - <span style="font-size: 14px"> S Zhou, L Pan, N Xiu,  and G  Li, A 0/1 constrained optimization solving sample average approximation for chance constrained programming, MOR, 2024. </span>

<!--
- <a style="font-size: 14px;color:#000000" href="https://jmlr.org/papers/v22/19-026.html" target="_blank"> S Zhou, N Xiu and H  Qi, Global and quadratic convergence of Newton hard-thresholding pursuit, *J Mach Learn Res*, 22:1−45, 2021.</a>
- <a style="font-size: 14px;color:#000000" href="https://www.sciencedirect.com/science/article/pii/S1063520322000458" target="_blank"> S Zhou, Gradient projection newton pursuit for sparsity constrained optimization, *Appl Comput Harmon Anal*, 61:75-100, 2022.</a> 
- <a style="font-size: 14px;color:#000000" href="http://www.yokohamapublishers.jp/online2/oppjo/vol13/p325.html" target="_blank"> L Pan, S Zhou, N Xiu, and H Qi, A convergent iterative hard thresholding for nonnegative sparsity optimization, *Pac J Optim*, 13:325-353, 2017.</a>  

---
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
 -->
