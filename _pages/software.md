---
layout: archive
title: ""  
permalink: /software/
author_profile: true
redirect_from:
  - /resume
---
<details>
  <summary><span style="color:grey"><b style="font-size:25px">0/1 Loss Optimization Solvers</b></span></summary>

  <a href="https://github.com/ShenglongZhou/GPSP">GPSP</a>, a matlab package solving the one-bit compressive sensing problems.  Source codes for <br>
  <a href="https://www.researchgate.net/publication/348371863">Computing One-bit Compressive Sensing via Double-Sparsity Constrained Optimization</a>.
  <br> 
  <br>
  <a href="https://github.com/Huajun-Wang/L01ADMM">L01ADMM</a>, a matlab package solving the support vector machine.  Source codes for <br>
  <a href="https://arxiv.org/abs/1912.07418">Support vector machine classifier via  $L_{0/1}$ soft-margin loss</a>.  
  
</details> 

<details>
  <summary><span style="color:grey"><b style="font-size:25px">Sparse  Optimization Solvers</b></span></summary>
 
  <a href="https://github.com/ShenglongZhou/NHTPver2">NHTP</a>, a matlab package solving the sparsity constrained optimization problems including compressed sensing, sparse logistic regression, sparse linear complementarity problems and so on.  Source codes for <br>
  <a href="https://arxiv.org/abs/1901.02763">Global and Quadratic Convergence of Newton Hard-Thresholding Pursuit</a>.
  
  <a href="https://github.com/ShenglongZhou/NL0R">NL0R</a>, a matlab package solving the $\ell_0$ regularized optimization problems including compressed sensing, sparse logistic regression, sparse linear complementarity problems.  Source codes for <br>
  <a href="https://arxiv.org/abs/2004.05132"Newton Method for $\ell_0$-Regularized Optimization</a>>.

  <a href="https://github.com/ShenglongZhou/IIHT">IIHT</a>, a matlab package solving the sparsity constrained optimization problems.  Source codes for <br>
  <a href="http://www.ybook.co.jp/online2/oppjo/vol13/p325.html">A Convergent Iterative Hard Thresholding for Sparsity and Nonnegativity Constrained Optimization</a>. 

  <a href="https://github.com/ShenglongZhou/MIRL1">MIRL1</a>, a matlab package solving the reweighted $\ell_1$ minimization.  Source codes for <br>
   <a href="https://doi.org/10.1093/imaiai/iaw002">A Null-space-based Weighted $\ell_1$ Minimisation Approach to Compressed Sensing</a>.

  <a href="https://github.com/ShenglongZhou/HTPCP">HTPCP</a>, a matlab package solving the sparse linear/nonlinear complementarity problems.  Source codes for <br>
  <a href="https://link.springer.com/article/10.1007/s11590-014-0834-7">A Half Thresholding Projection Algorithmfor Sparse Solutions of LCPs</a>. 

  <a href="https://github.com/ShenglongZhou/NSSVM">NSSVM</a>, a matlab package solving the sparse support vector machine.  Source codes for <br>
  <a href="https://arxiv.org/abs/2005.13771">Sparse SVM for Sufficient Data Reduction</a>. 

  <a href="https://github.com/ShenglongZhou/ADMM">ADMM</a>, a matlab package solving the  sparse and low-rank covariance matrix recovery problem.  Source codes for <br>
  <a href="https://link.springer.com/article/10.1007/s40305-014-0058-7">Sparse and Low-Rank Covariance Matrix Estimation</a>. 
 

#### Two general forms of sparse optimization 
 ---
  * Sparsity  constrained optimization:
\begin{eqnarray}
\label{SCO} \min_{x} && f(x), ~ {\rm s.t.}, ~ \Vert x \Vert_0\leq s 
\end{eqnarray}
 where $f: \mathbb{R}^{ n}\rightarrow  \mathbb{R}$, $s\ll n$ and $\Vert x \Vert_0$ is the so-called $\ell_0$ norm that counts the number of nonzero elements of $x$. 
 
 * $\ell_0$ regularized optimization:
\begin{eqnarray}
\label{L0RO} \min_{x} && f(x) +\lambda \Vert x \Vert_0 
\end{eqnarray}
 where $f: \mathbb{R}^{ n}\rightarrow  \mathbb{R}$ and $\lambda>0$. 
 
#### Applications of sparse optimization:  
---

 * Compressed sensing (<span style="color:orange">**CS**</span>):
\begin{eqnarray}
f(x) = (1/2) \Vert Ax-b \Vert^2
\end{eqnarray}
where $A\in\mathbb{R}^{m\times n}, b\in \mathbb{R}^{m}$. 

* Sparse logistic regression (<span style="color:orange">**SLR**</span>):
\begin{eqnarray}
f(x) =  \frac{1}{m}\sum_{i=1}^{m}\left\lbrace \ln(1+ e^{\langle a_i, x\rangle})-b_i\langle a_i, x\rangle\right \rbrace+\mu\Vert x\Vert_2^2  
\end{eqnarray}
where $a_i\in\mathbb{R}^{n}, b_i\in \lbrace 0,1\rbrace, i=1,2,\cdots,m$ and $\mu\geq0$.

* Sparse linear complementarity problem (<span style="color:orange">**SLCP**</span>):
\begin{eqnarray}
f(x) = \frac{1}{r}\sum_{i=1}^{m}\left\lbrace   (x_i)^r_{+}(M_ix+q_i)^r_{+}  +   (-x_i)^r_{+}   +  (-M_ix-q_i)^r_+ \right \rbrace 
\end{eqnarray}
where $M\in\mathbb{R}^{n\times n}, q\in \mathbb{R}^{n}, r\geq 2$, $M_i$ is the $i$th row of $M$ and $t_+:=\max \lbrace t,0\rbrace$. 
Note that  
\begin{eqnarray}
 f(x)=0~~ \Longleftrightarrow~~ x \geq 0,~ Mx+q\geq 0,~ \langle x , Mx+q \rangle=0 \nonumber
\end{eqnarray}

Applications solved by the aforementioned solvers are summarized in following table:

 <table border="2" width="0.5">
    <tr>
      <td style="width:5%" align="center"> </td>
      <td style="width:5%" align="center"><a  href='https://github.com/ShenglongZhou/NHTPver2'>NHTP</a></td>
      <td style="width:5%" align="center"><a  href='https://github.com/ShenglongZhou/NL0R'>NL0R</a></td>
      <td style="width:5%" align="center"><a  href='https://github.com/ShenglongZhou/IIHT'>IIHT</a></td>
      <td style="width:5%" align="center"><a  href='https://github.com/ShenglongZhou/MIRL1'>MIRL1</a></td>
      <td style="width:5%" align="center"><a  href='https://github.com/ShenglongZhou/HTPCP'>HTPCP</a></td>
    </tr>
     <tr>
    	  <td style="width:5%" align="left"><span style="color:orange">${\bf {\rm CS}}$</span></td>
        <td style="width:5%" align="center">$\surd$</td>
        <td style="width:5%" align="center">$\surd$</td>
        <td style="width:5%" align="center">$\surd$</td>
        <td style="width:5%" align="center">$\surd$</td>
        <td style="width:5%" align="center"> </td> 
    </tr>
      <tr>
    	  <td style="width:5%" align="left"><span style="color:orange">${\bf {\rm SLR}}$</span></td>
        <td style="width:5%" align="center">$\surd$</td>
        <td style="width:5%" align="center">$\surd$</td>
        <td style="width:5%" align="center">$\surd$</td>
        <td style="width:5%" align="center"> </td> 
        <td style="width:5%" align="center"> </td> 
    </tr>
      <tr>
    	  <td style="width:5%" align="left"><span style="color:orange">${\bf {\rm SLCP}}$</span></td>
        <td style="width:5%" align="center">$\surd$</td>
        <td style="width:5%" align="center">$\surd$</td>
        <td style="width:5%" align="center">$\surd$</td>
        <td style="width:5%" align="center"> </td>
        <td style="width:5%" align="center">$\surd$</td> 
    </tr>
    </table>
</details>

<details>
  <summary><span style="color:grey"><b style="font-size:25px">EDM Optimization Solvers</b></span></summary>
 
* [SQREDM](https://github.com/ShenglongZhou/SQREDM), a matlab package solving Euclidean distance matrix optimization problems for multidimensional scaling including sensor network localization, molecular conformation.   Source codes for <br>
  [A Fast Matrix Majorization-Projection Method for Penalized Stress Minimization with Box Constraints](https://ieeexplore.ieee.org/document/8399531).

* [PREEEDM](https://github.com/ShenglongZhou/PREEEDM), a matlab package solving Euclidean distance matrix optimization problems for multidimensional scaling including sensor network localization, molecular conformation.   Source codes for <br>
  [Robust Euclidean Embedding via EDM Optimization](https://doi.org/10.1007/s12532-019-00168-0).
</details>

 
<details>
  <summary><span style="color:grey"><b style="font-size:25px">Bilevel Optimization Toolbox </b></span></summary>
  

* [BOLIB](https://github.com/ShenglongZhou/BOLIB), the first version of the library providing 124 test examples. See [BOLIB: bilevel optimization library of test problems](https://arxiv.org/abs/1812.00230) for more information.

* [BOLIBver2](https://biopt.github.io/bolib/), the second version of  the library providing 173 test examples. See [BOLIB2019: bilevel optimization library of test problems version 2](https://www.researchgate.net/publication/338375731) for more information.

* [BiOpt Toolbox](https://biopt.github.io/),  bilevel optimization toolbox including 173 test examples from [BOLIBver2](https://biopt.github.io/bolib/), three solvers to solve bilevel optimization problems and several useful tools.
</details>




