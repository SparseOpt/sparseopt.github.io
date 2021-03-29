---
layout: archive
title: ""  
permalink: /software/
author_profile: true
redirect_from:
  - /resume
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


## <span style="color:#225675"><b style="font-size:20px">0/1 Loss Optimization Solvers</b></span>
---

  <font size=4> 
  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/ShenglongZhou/GPSP" target="_blank">GPSP</a> solves the one-bit compressive sensing problems. <br>
  <a style="color:#225675" href="https://www.researchgate.net/publication/348371863" target="_blank"><i>Computing One-bit Compressive Sensing via Double-Sparsity Constrained Optimization</i></a>. <br> <br>
  
  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/Huajun-Wang/L01ADMM" target="_blank">L01ADMM</a> solves the support vector machine. <br>
  <a style="color:#225675" href="https://arxiv.org/abs/1912.07418" target="_blank"><i>Support vector machine classifier via  L0/1 soft-margin loss</i></a>. 
      
  </font>

## <span style="color:#225675"><b style="font-size:20px">Sparse  Optimization Solvers</b></span>
---

  <font size=4> 
  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/ShenglongZhou/NHTPver2" target="_blank">NHTP</a> solves the sparsity constrained optimization including CS, LR, LCP and etc.   <br> 
  <a style="color:#225675" href="https://arxiv.org/abs/1901.02763" target="_blank"><i>Global and Quadratic Convergence of Newton Hard-Thresholding Pursuit</i></a>. <br> <br>
    
  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/ShenglongZhou/NL0R" target="_blank">NL0R</a> solves the L0 regularized optimization problems including CS, LR and LCP.   <br> 
  <a style="color:#225675" href="https://arxiv.org/abs/2004.05132" target="_blank"><i>Newton Method for L0-Regularized Optimization</i></a>.<br>  <br> 
  
  <details>
  <summary><span style="color:#007D98"><b style="font-size:15px">Click for more solvers</b></span></summary>
  <br> 

  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/ShenglongZhou/IIHT" target="_blank">IIHT</a> solves the sparsity constrained optimization including CS, LR, LCP and etc.   <br>
  <a style="color:#225675" href="http://www.ybook.co.jp/online2/oppjo/vol13/p325.html" target="_blank"><i>A Convergent Iterative Hard Thresholding for Sparsity and Nonnegativity Constrained Optimization</i></a>. <br><br>
 
  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/ShenglongZhou/MIRL1" target="_blank">MIRL1</a> solves the reweighted L1 minimization.    <br>
  <a style="color:#225675" href="https://doi.org/10.1093/imaiai/iaw002" target="_blank"><i>A Null-space-based Weighted L1 Minimisation Approach to Compressed Sensing</i></a>.<br><br>
 
  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/ShenglongZhou/HTPCP" target="_blank">HTPCP</a> solves the sparse linear/nonlinear complementarity problems.   <br>
  <a style="color:#225675" href="https://link.springer.com/article/10.1007/s11590-014-0834-7" target="_blank"><i>A Half Thresholding Projection Algorithmfor Sparse Solutions of LCPs</i></a>. <br><br>
 
  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/ShenglongZhou/NSSVM" target="_blank">NSSVM</a> solves the sparse support vector machine.  Source codes for <br>
  <a style="color:#225675" href="https://arxiv.org/abs/2005.13771" target="_blank"><i>Sparse SVM for Sufficient Data Reduction</i></a>. <br><br>
 
  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/ShenglongZhou/ADMM" target="_blank">ADMM</a> solves the sparse and low-rank covariance matrix recovery problem.   <br>
  <a style="color:#225675" href="https://link.springer.com/article/10.1007/s40305-014-0058-7" target="_blank"><i>Sparse and Low-Rank Covariance Matrix Estimation</i></a>. <br><br>
 

  <b> Two general forms of sparse optimization: </b> 
  
   \begin{eqnarray*}
   \begin{array}{lll}
   \text{Sparsity constrained optimization:}~&~\min_{x}&~f(x), ~ {\rm s.t.}, ~ \Vert x \Vert_0\leq s \\
   \text{L0 regularized optimization:} &~\min_{x}&~f(x) +\lambda \Vert x \Vert_0,
   \end{array}
   \end{eqnarray*}
   where $f: \mathbb{R}^{ n}\rightarrow  \mathbb{R}$, $s\ll n, \lambda>0$ and $\Vert x \Vert_0$ is the so-called $\ell_0$ norm that counts the number of nonzero elements of $x$.  <br><br>
 
 
  <!---### <b> Applications of sparse optimization </b>  <br><br>
  * Compressed sensing (<span style="color:orange"><b>CS</b></span>):
  \begin{eqnarray}
  f(x) = (1/2) \Vert Ax-b \Vert^2
  \end{eqnarray}
  where $A\in\mathbb{R}^{m\times n}, b\in \mathbb{R}^{m}$. <br><br> 
  * Sparse logistic regression (<span style="color:orange"><b>SLR</b></span>):
  \begin{eqnarray}
  f(x) =  \frac{1}{m}\sum_{i=1}^{m}\left\lbrace \ln(1+ e^{\langle a_i, x\rangle})-b_i\langle a_i, x\rangle\right \rbrace+\mu\Vert x\Vert_2^2  
  \end{eqnarray}
  where $a_i\in\mathbb{R}^{n}, b_i\in \lbrace 0,1\rbrace, i=1,2,\cdots,m$ and $\mu\geq0$.<br><br>
  * Sparse linear complementarity problem (<span style="color:orange"><b>SLCP</b></span>):
  \begin{eqnarray}
  f(x) = \frac{1}{r}\sum_{i=1}^{m}\left\lbrace   (x_i)^r_{+}(M_ix+q_i)^r_{+}  +   (-x_i)^r_{+}   +  (-M_ix-q_i)^r_+ \right \rbrace 
  \end{eqnarray}
  where $M\in\mathbb{R}^{n\times n}, q\in \mathbb{R}^{n}, r\geq 2$, $M_i$ is the $i$th row of $M$ and $t_+:=\max \lbrace t,0\rbrace$. 
  Note that  
  \begin{eqnarray}
   f(x)=0~~ \Longleftrightarrow~~ x \geq 0,~ Mx+q\geq 0,~ \langle x , Mx+q \rangle=0 \nonumber
  \end{eqnarray}
  <br>
  --->

   Applications solved by the aforementioned solvers are summarized in following table:<br>

   <table border="2" width="0.5">
      <tr>
        <td style="width:8%" align="center"> </td>
        <td style="width:5%" align="center"><a style="font-size: 16px; font-weight: bold; color:#007D98"  href="https://github.com/ShenglongZhou/NHTPver2" target="_blank">NHTP</a></td>
        <td style="width:5%" align="center"><a style="font-size: 16px; font-weight: bold; color:#007D98"  href="https://github.com/ShenglongZhou/NL0R" target="_blank">NL0R</a></td>
        <td style="width:5%" align="center"><a style="font-size: 16px; font-weight: bold; color:#007D98"  href="https://github.com/ShenglongZhou/IIHT" target="_blank">IIHT</a></td>
        <td style="width:5%" align="center"><a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/ShenglongZhou/MIRL1" target="_blank">MIRL1</a></td>
        <td style="width:5%" align="center"><a style="font-size: 16px; font-weight: bold; color:#007D98"  href="https://github.com/ShenglongZhou/HTPCP" target="_blank">HTPCP</a></td>
      </tr>
       <tr>
          <td style="width:8%" align="left"><b>Compressed sensing (<span style="color:#007D98">CS</span>)</b></td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center"> </td> 
      </tr>
        <tr>
          <td style="width:8%" align="left"><b>Logistic regression (<span style="color:#007D98">LR</span>)</b></td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center"> </td> 
          <td style="width:5%" align="center"> </td> 
      </tr>
        <tr>
          <td style="width:8%" align="left"><b>Linear complementarity problem (<span style="color:#007D98">LCP</span>)</b></td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center"> </td>
          <td style="width:5%" align="center">$\surd$</td> 
      </tr>
      </table>
  </details> 
  </font>


## <span style="color:#225675"><b style="font-size:20px">EDM Optimization Solvers</b></span>
---
  
  <font size=4> 

  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/ShenglongZhou/SQREDM" target="_blank">SQREDM</a> solves applications from multidimensional scaling,
  e.g.,  sensor network localization and molecular conformation.    <br>
  <a style="color:#225675" href="https://ieeexplore.ieee.org/document/8399531" target="_blank"><i>A Fast Matrix Majorization-Projection Method for Penalized Stress Minimization with Box Constraints</i></a>.<br><br>
 
  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://github.com/ShenglongZhou/PREEEDM" target="_blank">PREEEDM</a> solves applications from multidimensional scaling,
  e.g.,  sensor network localization and molecular conformation.    <br>
  <a style="color:#225675" href="https://doi.org/10.1007/s12532-019-00168-0" target="_blank"><i>Robust Euclidean Embedding via EDM Optimization</i></a>. 
  </font>



## <span style="color:#225675"><b style="font-size:20px">Bilevel Optimization Toolbox </b></span>
---

  <font size=4>
 
  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://biopt.github.io/bolib/" target="_blank">BOLIBver2</a>, the second version of  the library providing 173 test examples. <br>
  <a style="color:#225675" href="https://www.researchgate.net/publication/338375731" target="_blank"><i>BOLIB2019: bilevel optimization library of test problems version 2</i></a>.<br><br>
 
  <a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://biopt.github.io/" target="_blank">BiOpt Toolbox</a>,  bilevel optimization toolbox including <a style="color:#007D98" href="https://biopt.github.io/bolib/" target="_blank">BOLIBver2</a>, 
  <a style="color:#225675" href="https://biopt.github.io/solvers/" target="_blank">three solvers</a> and several useful tools. 
  
  </font>
