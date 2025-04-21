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
where  $f:\mathbb{R}^{n}\rightarrow \mathbb{R}$ is a twice continuously differentiable function,  $\mathbf{A}\in\mathbb{R}^{m\times n}$ is a matrix, $\mathbf{b}\in\mathbb{R}^{m}$ is a vector, $\lambda$ is a given scalar,  $(\mathbf{z})_+=(\max\{0,z_1\},\ldots,\max\{0,z_m\})^\top$, and  $\|\mathbf{z}\|_0$ is the so-called $\ell_0$ norm that counts the number of nonzero entries in $\mathbf{z}$. The regularization term is related to the step function (or 0/1 loss function) defined by $\ell_{0/1}(t)=1$ if $t>0$ and $\ell_{0/1}(t)=0$ otherwise. Therefore,
  \begin{equation}\|(\mathbf{z})_+\|_0=\sum_{i=1}^m \ell_{0/1}\left(z_{i}\right)\nonumber\end{equation}
</div>
 
<!-- ## <span style="color:#8C8C8C"> The solver and its demonstration </span> -->

---
<div style="text-align:justify;"> 
The solver can be download here - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href="https://github.com/ShenglongZhou/NM01" target="_blank">NM01</a>,
which was developed from the following paper:
</div>

> <span style="font-size: 14px"> S Zhou, L Pan, N Xiu,  and H Qi, Quadratic convergence of smoothing Newton's method for 0/1 loss optimization, SIAM J Optim, 31:3184–3211, 2021. </span>

<!--
- <a style="font-size: 14px;color:#000000" href="https://jmlr.org/papers/v22/19-026.html" target="_blank"> S Zhou, N Xiu and H  Qi, Global and quadratic convergence of Newton hard-thresholding pursuit, *J Mach Learn Res*, 22:1−45, 2021.</a>
- <a style="font-size: 14px;color:#000000" href="https://www.sciencedirect.com/science/article/pii/S1063520322000458" target="_blank"> S Zhou, Gradient projection newton pursuit for sparsity constrained optimization, *Appl Comput Harmon Anal*, 61:75-100, 2022.</a> 
- <a style="font-size: 14px;color:#000000" href="http://www.yokohamapublishers.jp/online2/oppjo/vol13/p325.html" target="_blank"> L Pan, S Zhou, N Xiu, and H Qi, A convergent iterative hard thresholding for nonnegative sparsity optimization, *Pac J Optim*, 13:325-353, 2017.</a> -->

---
<div style="text-align:justify;">  
Note that <b style="font-size:14px;color:#777777">NM01</b> is a second-order method, which require the gradient and Hessian of $f$. Below is a demonstration of how to define the gradient and Hessian for the solver.
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
