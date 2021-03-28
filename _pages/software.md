---
layout: archive
title: ""  
permalink: /software/
author_profile: true
redirect_from:
  - /resume
---
 
 
## <span style="color:grey"><b style="font-size:20px">0/1 Loss Optimization Solvers</b></span>
---

  <font size=3> 
  <a href="https://github.com/ShenglongZhou/GPSP">GPSP</a> solves the one-bit compressive sensing problems. <br>
  <a href="https://www.researchgate.net/publication/348371863">Computing One-bit Compressive Sensing via Double-Sparsity Constrained Optimization</a>. <br> <br>
  
  <a href="https://github.com/Huajun-Wang/L01ADMM">L01ADMM</a> solves the support vector machine. <br>
  <a href="https://arxiv.org/abs/1912.07418">Support vector machine classifier via  $L_{0/1}$ soft-margin loss</a>. <br> <br>
      
  </font>



## <span style="color:grey"><b style="font-size:20px">Sparse  Optimization Solvers</b></span>
---

  <font size=3> 
  <a href="https://github.com/ShenglongZhou/NHTPver2">NHTP</a> solves the sparsity constrained optimization, 
  e.g., compressed sensing, logistic regression, linear complementarity problems and etc.   <br> 
  <a href="https://arxiv.org/abs/1901.02763">Global and Quadratic Convergence of Newton Hard-Thresholding Pursuit</a>. <br> <br>
    
  <a href="https://github.com/ShenglongZhou/NL0R">NL0R</a> solves the $\ell_0$ regularized optimization problems,
  e.g., compressed sensing, logistic regression,  linear complementarity problems.   <br> 
  <a href="https://arxiv.org/abs/2004.05132">Newton Method for $\ell_0$-Regularized Optimization</a>.<br>  <br> 
  
  <details>
  <summary><span style="color:grey"><b style="font-size:10px">Click for more solvers</b></span></summary>
  <br> 

  <a href="https://github.com/ShenglongZhou/IIHT">IIHT</a>, a matlab package solving the sparsity constrained optimization problems.   <br>
  <a href="http://www.ybook.co.jp/online2/oppjo/vol13/p325.html">A Convergent Iterative Hard Thresholding for Sparsity and Nonnegativity Constrained Optimization</a>. <br><br>
 
  <a href="https://github.com/ShenglongZhou/MIRL1">MIRL1</a>, a matlab package solving the reweighted $\ell_1$ minimization.    <br>
  <a href="https://doi.org/10.1093/imaiai/iaw002">A Null-space-based Weighted $\ell_1$ Minimisation Approach to Compressed Sensing</a>.<br><br>
 
  <a href="https://github.com/ShenglongZhou/HTPCP">HTPCP</a>, a matlab package solving the sparse linear/nonlinear complementarity problems.   <br>
  <a href="https://link.springer.com/article/10.1007/s11590-014-0834-7">A Half Thresholding Projection Algorithmfor Sparse Solutions of LCPs</a>. <br><br>
 
  <a href="https://github.com/ShenglongZhou/NSSVM">NSSVM</a>, a matlab package solving the sparse support vector machine.  Source codes for <br>
  <a href="https://arxiv.org/abs/2005.13771">Sparse SVM for Sufficient Data Reduction</a>. <br><br>
 
  <a href="https://github.com/ShenglongZhou/ADMM">ADMM</a>, a matlab package solving the  sparse and low-rank covariance matrix recovery problem.   <br>
  <a href="https://link.springer.com/article/10.1007/s40305-014-0058-7">Sparse and Low-Rank Covariance Matrix Estimation</a>. <br><br>
 

  <b> Two general forms of sparse optimization </b> <br><br>

  Sparsity  constrained optimization:
  \begin{eqnarray}
  \label{SCO} \min_{x} && f(x), ~ {\rm s.t.}, ~ \Vert x \Vert_0\leq s 
  \end{eqnarray}
   where $f: \mathbb{R}^{ n}\rightarrow  \mathbb{R}$, $s\ll n$ and $\Vert x \Vert_0$ is the so-called $\ell_0$ norm that counts the number of nonzero elements of $x$. 

  <br> <br>
  $\ell_0$ regularized optimization:
  \begin{eqnarray}
  \label{L0RO} \min_{x} && f(x) +\lambda \Vert x \Vert_0 
  \end{eqnarray}
   where $f: \mathbb{R}^{ n}\rightarrow  \mathbb{R}$ and $\lambda>0$.  <br><br>
 
 
  <b> Applications of sparse optimization </b>  <br><br>
 
  Compressed sensing (<span style="color:orange"><b>CS</b></span>):
  \begin{eqnarray}
  f(x) = (1/2) \Vert Ax-b \Vert^2
  \end{eqnarray}
  where $A\in\mathbb{R}^{m\times n}, b\in \mathbb{R}^{m}$. <br><br>

 
  Sparse logistic regression (<span style="color:orange"><b>SLR</b></span>):
  \begin{eqnarray}
  f(x) =  \frac{1}{m}\sum_{i=1}^{m}\left\lbrace \ln(1+ e^{\langle a_i, x\rangle})-b_i\langle a_i, x\rangle\right \rbrace+\mu\Vert x\Vert_2^2  
  \end{eqnarray}
  where $a_i\in\mathbb{R}^{n}, b_i\in \lbrace 0,1\rbrace, i=1,2,\cdots,m$ and $\mu\geq0$.<br><br>


  Sparse linear complementarity problem (<span style="color:orange"><b>SLCP</b></span>):
  \begin{eqnarray}
  f(x) = \frac{1}{r}\sum_{i=1}^{m}\left\lbrace   (x_i)^r_{+}(M_ix+q_i)^r_{+}  +   (-x_i)^r_{+}   +  (-M_ix-q_i)^r_+ \right \rbrace 
  \end{eqnarray}
  where $M\in\mathbb{R}^{n\times n}, q\in \mathbb{R}^{n}, r\geq 2$, $M_i$ is the $i$th row of $M$ and $t_+:=\max \lbrace t,0\rbrace$. 
  Note that  
  \begin{eqnarray}
   f(x)=0~~ \Longleftrightarrow~~ x \geq 0,~ Mx+q\geq 0,~ \langle x , Mx+q \rangle=0 \nonumber
  \end{eqnarray}
  <br> 

  Applications solved by the aforementioned solvers are summarized in following table:<br>

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
          <td style="width:5%" align="left"><span style="color:orange"><b>CS</b></span></td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center"> </td> 
      </tr>
        <tr>
          <td style="width:5%" align="left"><span style="color:orange"><b>SLR</b></span></td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center"> </td> 
          <td style="width:5%" align="center"> </td> 
      </tr>
        <tr>
          <td style="width:5%" align="left"><span style="color:orange"><b>SLCP</b></span></td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center"> </td>
          <td style="width:5%" align="center">$\surd$</td> 
      </tr>
      </table>
  </details> 
  </font>


## <span style="color:grey"><b style="font-size:25px">EDM Optimization Solvers</b></span>
---
  
  <font size=3> 

  <a href="https://github.com/ShenglongZhou/SQREDM">SQREDM</a> solves EDM optimization for multidimensional scaling,
  e.g.,  sensor network localization and molecular conformation.    <br>
  <a href="https://ieeexplore.ieee.org/document/8399531">A Fast Matrix Majorization-Projection Method for Penalized Stress Minimization with Box Constraints</a>.<br><br>
 
  <a href="https://github.com/ShenglongZhou/PREEEDM">PREEEDM</a> solves EDM optimization for multidimensional scaling ,
  e.g.,  sensor network localization and molecular conformation.    <br>
  <a href="https://doi.org/10.1007/s12532-019-00168-0">Robust Euclidean Embedding via EDM Optimization</a>. 
  </font>



## <span style="color:grey"><b style="font-size:25px">Bilevel Optimization Toolbox </b></span>
---

  <font size=3>
 
  <a href="https://github.com/ShenglongZhou/BOLIB">BOLIB</a>, the first version of the library providing 124 test examples.  <br>
  <a href="https://arxiv.org/abs/1812.00230">BOLIB: bilevel optimization library of test problems</a>.<br><br>
 
  <a href="https://biopt.github.io/bolib/">BOLIBver2</a>, the second version of  the library providing 173 test examples. <br>
  <a href="https://www.researchgate.net/publication/338375731">BOLIB2019: bilevel optimization library of test problems version 2</a>.<br><br>
 
  <a href="https://biopt.github.io/">BiOpt Toolbox</a>,  bilevel optimization toolbox including <a href="https://biopt.github.io/bolib/">BOLIBver2</a>,  
  <a href="https://biopt.github.io/solvers/">three solvers</a>and several useful tools. 
  </font>


