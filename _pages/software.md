---
layout: archive
title: ""  
permalink: /software/
author_profile: true
redirect_from:
  - /resume
---

<span style="color:grey">Bilevel optimization toolbox</span> 
---

* [BOLIB](https://github.com/ShenglongZhou/BOLIB), the first version of the library providing 124 test examples. See [BOLIB: bilevel optimization library of test problems](https://arxiv.org/abs/1812.00230) for more information.

* [BOLIBver2](https://biopt.github.io/bolib/), the second version of  the library providing 173 test examples. See [BOLIB2019: bilevel optimization library of test problems version 2](https://www.researchgate.net/publication/338375731) for more information.

* [BiOpt Toolbox](https://biopt.github.io/),  bilevel optimization toolbox including 173 test examples from [BOLIBver2](https://biopt.github.io/bolib/), three solvers to solve bilevel optimization problems and several useful tools.

<span style="color:grey">Euclidean distance matrix optimization solvers</span> 
---

* [SQREDM](https://github.com/ShenglongZhou/SQREDM), a matlab package solving Euclidean distance matrix optimization problems for multidimensional scaling including sensor network localization, molecular conformation.   Source codes for <br>
  [A Fast Matrix Majorization-Projection Method for Penalized Stress Minimization with Box Constraints](https://ieeexplore.ieee.org/document/8399531).

* [PREEEDM](https://github.com/ShenglongZhou/PREEEDM), a matlab package solving Euclidean distance matrix optimization problems for multidimensional scaling including sensor network localization, molecular conformation.   Source codes for <br>
  [Robust Euclidean Embedding via EDM Optimization](https://doi.org/10.1007/s12532-019-00168-0).



<span style="color:grey">Sparse optimization solvers</span> 
---

* [NHTP](https://github.com/ShenglongZhou/NHTP) or [NHTPver2](https://github.com/ShenglongZhou/NHTPver2), a matlab package solving the sparsity constrained optimization problems including compressed sensing, sparse logistic regression, sparse linear complementarity problems and so on.  Source codes for <br>
  [Global and Quadratic Convergence of Newton Hard-Thresholding Pursuit](https://arxiv.org/abs/1901.02763).

* [NL0R](https://github.com/ShenglongZhou/NL0R), a matlab package solving the $\ell_0$ regularized optimization problems including compressed sensing, sparse logistic regression, sparse linear complementarity problems.  Source codes for <br>
  [Newton Method for $\ell_0$-Regularized Optimization](https://arxiv.org/abs/2004.05132).

* [IIHT](https://github.com/ShenglongZhou/IIHT), a matlab package solving the sparsity constrained optimization problems.  Source codes for <br>
[A Convergent Iterative Hard Thresholding for Sparsity and Nonnegativity Constrained Optimization](http://www.ybook.co.jp/online2/oppjo/vol13/p325.html). 

* [MIRL1](https://github.com/ShenglongZhou/MIRL1), a matlab package solving the reweighted $\ell_1$ minimization.  Source codes for <br>
  [A Null-space-based Weighted $\ell_1$ Minimisation Approach to Compressed Sensing](https://doi.org/10.1093/imaiai/iaw002).
 
* [HTPCP](https://github.com/ShenglongZhou/HTPCP), a matlab package solving the sparse linear/nonlinear complementarity problems.  Source codes for <br>
  [A Half Thresholding Projection Algorithmfor Sparse Solutions of LCPs](https://link.springer.com/article/10.1007/s11590-014-0834-7). 
  
* [NSSVM](https://github.com/ShenglongZhou/NSSVM), a matlab package solving the sparse support vector machine.  Source codes for <br>
  [Sparse SVM for Sufficient Data Reduction](https://arxiv.org/abs/2005.13771). 
  
* [ADMM](https://github.com/ShenglongZhou/ADMM), a matlab package solving the  sparse and low-rank covariance matrix recovery problem.  Source codes for <br>
  [Sparse and Low-Rank Covariance Matrix Estimation](https://link.springer.com/article/10.1007/s40305-014-0058-7). 
 
### Two general forms of Sparse optimization 
 
  * Sparsity  constrained optimization (<span style="color:orange">**SCO**</span>):
\begin{eqnarray}
\label{SCO} \min_{x} && f(x), ~ {\rm s.t.}, ~ \Vert x \Vert_0\leq s.
\end{eqnarray}
 where $f: \mathbb{R}^{ n}\rightarrow  \mathbb{R}$, $s\ll n$ and $\Vert x \Vert_0$ is the so-called $\ell_0$ norm that counts the number of nonzero elements of $x$. 
 
 *  $\ell_0$ regularized optimization (<span style="color:orange">**L0RO**</span>):
\begin{eqnarray}
\label{L0RO} \min_{x} && f(x) +\lambda \Vert x \Vert_0\leq s.
\end{eqnarray}
 where $f: \mathbb{R}^{ n}\rightarrow  \mathbb{R}$, $\lambda>0$. 
 
### Applications of Sparse optimization  
It has a wide range of applications including:
 * Compressed sensing (<span style="color:orange">**CS**</span>):
\begin{eqnarray}
\label{CSS}\min_{x} && \Vert Ax-b \Vert^2, ~ {\rm s.t.}, ~ \Vert x \Vert_0\leq s,  \\\\\\
\label{CSL}\min_{x} && \Vert Wx \Vert_1, ~~~~~~~ \mbox{s.t.},~ Ax=b, 
\end{eqnarray}
where $A\in\mathbb{R}^{m\times n}, b\in \mathbb{R}^{m}, W\in\mathbb{R}^{n\times n}$ is a diagonal matrix with $W_{ii}>0$ and $\Vert \cdot\Vert$ and $\Vert \cdot\Vert_1$ are the $\ell_2$ and $\ell_1$ norm. 

* Sparse logistic regression (<span style="color:orange">**SLR**</span>):
\begin{eqnarray}
\label{SLR} \min_{x}~  \frac{1}{m}\sum_{i=1}^{m}\left\lbrace \ln(1+ e^{\langle a_i, x\rangle})-b_i\langle a_i, x\rangle\right \rbrace+\mu\Vert x\Vert_2^2 , ~ ~ {\rm s.t.},~ \Vert x\Vert_0\leq s.
\end{eqnarray}
where $a_i\in\mathbb{R}^{n}, b_i\in \lbrace 0,1\rbrace, i=1,2,\cdots,m$ and $\mu\geq0$.

* Sparse linear complementarity problem (<span style="color:orange">**SLCP**</span>):
\begin{eqnarray}
\label{SLCP} x \geq 0,~ Mx+q\geq 0,~ \langle x , Mx+q \rangle=0, ~ \Vert x\Vert_0\leq s.
\end{eqnarray}
where $M\in\mathbb{R}^{n\times n}$ and $q\in \mathbb{R}^{n}$. 

Applications solved by the above mentioned solvers are summarized in following table:

 <table border="2" width="0.5">
    <tr>
      <td style="width:10%" align="center"> </td>
      <td style="width:10%" align="center"><a  href='https://github.com/ShenglongZhou/NHTP'>NHTP</a></td>
      <td style="width:10%" align="center"><a  href='https://github.com/ShenglongZhou/NHTPver2'>NHTPver2</a></td>
      <td style="width:10%" align="center"><a  href='https://github.com/ShenglongZhou/IIHT'>IIHT</a></td>
      <td style="width:10%" align="center"><a  href='https://github.com/ShenglongZhou/MIRL1'>MIRL1</a></td>
      <td style="width:10%" align="center"><a  href='https://github.com/ShenglongZhou/HTPCP'>HTPCP</a></td>
    </tr>
    <tr>
    	  <td style="width:10%" align="left"><span style="color:orange">${\bf {\rm SCO}}$</span></td>
        <td style="width:10%" align="center">(\ref{SCO}) </td>
        <td style="width:10%" align="center">(\ref{SCO}) </td>
        <td style="width:10%" align="center">(\ref{SCO}) </td>
        <td style="width:10%" align="center"> </td>
        <td style="width:10%" align="center"> </td> 
    </tr>
     <tr>
    	  <td style="width:10%" align="left"><span style="color:orange">${\bf {\rm CS}}$</span></td>
        <td style="width:10%" align="center">(\ref{CSS}) </td>
        <td style="width:10%" align="center">(\ref{CSS}) </td>
        <td style="width:10%" align="center">(\ref{CSS}) </td>
        <td style="width:10%" align="center">(\ref{CSL})</td>
        <td style="width:10%" align="center"> </td> 
    </tr>
      <tr>
    	  <td style="width:10%" align="left"><span style="color:orange">${\bf {\rm SLR}}$</span></td>
        <td style="width:10%" align="center">(\ref{SLR}) </td>
        <td style="width:10%" align="center">(\ref{SLR}) </td>
        <td style="width:10%" align="center">(\ref{SLR}) </td>
        <td style="width:10%" align="center"> </td>
        <td style="width:10%" align="center"> </td> 
    </tr>
      <tr>
    	  <td style="width:10%" align="left"><span style="color:orange">${\bf {\rm SLCP}}$</span></td>
        <td style="width:10%" align="center"> </td>
        <td style="width:10%" align="center">(\ref{SLCP})</td>
        <td style="width:10%" align="center">(\ref{SLCP})</td>
        <td style="width:10%" align="center"> </td>
        <td style="width:10%" align="center">(\ref{SLCP})</td> 
    </tr>
    </table>
 
 



