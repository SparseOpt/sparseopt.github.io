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

 

## <span style="color:#8C8C8C">Sparsity-regularized optimization</span> 
---
<p style="line-height: 2;"></p>

\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n}} ~~  f(\mathbf{x}) + \lambda \parallel\mathbf{x}\parallel_q^q \tag{SRO}
\end{equation}

<div style="text-align:justify;">
where  $f:\mathbb{R}^{n}\rightarrow \mathbb{R}$ is a continuously or twice continuously differentiable function, $\lambda>0$ is a given scalar, and $\|\mathbf{x}\|_q^q=\sum_i |x_i|^q$ with $q\in[0,1)$ denotes the $\ell_q$-norm. In particular, when $q=0$,  $\|\mathbf{x}\|_0=\|\mathbf{x}\|_0^0$ is the so-called $\ell_0$-norm that counts the number of nonzero entries in $\mathbf{x}$.
</div>
 
<!-- ## <span style="color:#8C8C8C">The solver and its demonstration</span> -->
---
<div style="text-align:justify;">
The package can be download here - <a style="font-size: 16px; font-weight: bold; color:#006DB0" href="https://github.com/ShenglongZhou/CSpack" target="_blank">SROpack</a>,
which provides 2 solvers from the following papers:
</div>

> <b style="font-size:14px;color:#777777">PSNP</b> - <span style="font-size: 14px"> S Zhou, X Xiu, Y Wang, and D Peng, Revisiting Lq ($0 \leq q < 1$) norm regularized optimization, arXiv:2306.14394, 2023. </span>
<br><b style="font-size:14px;color:#777777">NL0R</b> - <span style="font-size: 14px"> S Zhou, L Pan, and N Xiu, Newton method for l0 regularized optimization, Numer Algorithms, 88:1541–1570, 2021. </span>

<!--
- <a style="font-size:14px; color:#000000" href="https://link.springer.com/article/10.1007/s11075-021-01085-x" target="_blank"> S Zhou, L Pan, and N Xiu, Newton method for L0-regularized optimization, *Numerical Algorithm*, 88:1541–1570, 2021 .</a>
- <a style="font-size:14px; color:#000000" href="https://arxiv.org/abs/2306.14394" target="_blank"> S Zhou, X Xiu, Y Wang, and D Peng, Revisiting Lq (0 <= q < 1) norm regularized optimization, *arXiv:2306.14394*, 2023.</a> 
-->

---
<div style="text-align:justify;">
Both solvers are second-order methods, which require both the gradient and Hessian of $f$. Below is a demonstration of how to define the gradient and Hessian for <b style="font-size:14px;color:#777777">NL0R</b> and <b style="font-size:14px;color:#777777">PSNP</b>.
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
