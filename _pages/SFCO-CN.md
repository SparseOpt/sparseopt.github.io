---
layout: archive
title: ""   
permalink: /SFCO-CN/
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

##  <span style="color:#8C8C8C"> Step function-constrained optimization</span> 
---

<p style="line-height: 1;"></p>
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{K}} ~~  f(\mathbf{x}),~~~~ \mbox{s.t.}~~ \parallel\mathbf{G}(\mathbf{x})\parallel_0^+\leq s,~\mathbf{x}\in \Omega  \tag{SFCO}
\end{equation}

<div style="text-align:justify;">
其中， $\mathbf{G}(\mathbf{x})\in\mathbb{R}^{M \times N}$ 是一个 $M \times N$ 维矩阵，每一个元素为 $G_{ij}(\mathbf{x})$，$1\leq i\leq M$，$1\leq j \leq N$，函数 $f:\mathbb{R}^{K}\rightarrow \mathbb{R}$ 和 $G_{ij}:\mathbb{R}^{K}\rightarrow \mathbb{R}$ 连续可微，最好是二次连续可微，$\Omega\subseteq\mathbb{R}^{K}$ 是一个闭凸集，$s$ 是一个远小于 $n$ 的正整数。 对于矩阵 $\mathbf{Z}\in\mathbb{R}^{M \times N}$，$\|\mathbf{Z}\|_0^+$ 计算该矩阵中含有正元素的列的个数，即
  \begin{equation}\|\mathbf{Z}\|_0^+= \mathrm{step}\Big(\max_{i=1,\ldots,M} Z_{i1}\Big)+\cdots+\mathrm{step}\Big(\max_{i=1,\ldots,M} Z_{iN}\Big)\nonumber\end{equation}
这里，$\mathrm{step}(t)$ 是阶梯函数，又称作 0/1 损失函数，定义为：如果 $t>0$，则 $\mathrm{step}(t)=1$；如果 $t\leq0$，则 $\mathrm{step}(t)=0$。特别地，当 $M=1$，矩阵 $\mathbf{Z}$ 退化成向量 $\mathbf{z}\in\mathbb{R}^{N}$，如果记 $\mathbf{z}_+=(\max\{0,z_1\},\ldots,\max\{0,z_N\})^\top$ 以及 $\parallel\mathbf{z}\parallel_0$ 表示 $\mathbf{z}$ 的零范数，也即 $\mathbf{z}$ 中非零元个数，则有  
  \begin{equation*}\|\mathbf{z}\|_0^+= \mathrm{step}(z_1)+\cdots+\mathrm{step}(z_N)=\|\mathbf{z}_+\|_0\end{equation*}
</div>
 
<!-- ## <span style="color:#8C8C8C"> The solver and its demonstration </span> -->

---
<div style="text-align:justify;"> 
程序包 - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href=" " target="_blank">SNSCO</a>（点击可直接下载）提供了 1 个求解器，其核心算法来自以下文章：
</div>

> <span style="font-size: 14px"> S Zhou, L Pan, N Xiu,  and G  Li, A 0/1 constrained optimization solving sample average approximation for chance constrained programming, MOR, 2024. </span>

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
