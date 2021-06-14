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




## <b style="font-size:20px">Sparse  Optimization Solvers</b>
---

  <font size=4> 
  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/CSpack" target="_blank">CSpack</a> offers four solvers to solve the compressive sensing problems.
  <p style="line-height: 1;"></p>
  
  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/NHTPver2" target="_blank">NHTP</a> solves the sparsity constrained optimization including CS, LR, LCP and etc. <a style="font-size: 16px; font-weight: bold; color:#8cd2d5"  href="https://www.jmlr.org/papers/volume22/19-026/19-026.pdf" target="_blank">Article</a> 
  <p style="line-height: 1;"></p>
    
  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/NL0R" target="_blank">NL0R</a> solves the L0 regularized optimization problems including CS, LR and LCP. <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://arxiv.org/abs/2004.05132" target="_blank">Article</a> 
  <p style="line-height: 1;"></p>

   <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/NSSVM" target="_blank">NSSVM</a> solves the sparse support vector machine.  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://arxiv.org/abs/2005.13771" target="_blank">Article</a> 
  <p style="line-height: 1;"></p>
  
  <!--- <details>
  <summary><span style="color:#8cd2d5"><b style="font-size:16px">Click for more solvers</b></span></summary>
  <br> --->

  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/IIHT" target="_blank">IIHT</a> solves the sparsity constrained optimization including CS, LR, LCP and etc. <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="http://www.ybook.co.jp/online2/oppjo/vol13/p325.html" target="_blank">Article</a>  
  <p style="line-height: 1;"></p>
 
  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/MIRL1" target="_blank">MIRL1</a> solves the reweighted L1 minimization.    <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://doi.org/10.1093/imaiai/iaw002" target="_blank">Article</a> 
  <p style="line-height: 1;"></p>
 
  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/HTPCP" target="_blank">HTPCP</a> solves the sparse linear/nonlinear complementarity problems. <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://link.springer.com/article/10.1007/s11590-014-0834-7" target="_blank">Article</a>  
  <p style="line-height: 1;"></p>

  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/ADMM" target="_blank">ADMM</a> solves the sparse and low-rank covariance matrix recovery problem.  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://link.springer.com/article/10.1007/s40305-014-0058-7" target="_blank">Article</a> <br><br>
 
   <!---
  <b> Two general forms of sparse optimization: </b> 
  
   \begin{eqnarray*}
   \begin{array}{lll}
   \text{Sparsity constrained optimization:}~&~\min_{x}&~f(x), ~ {\rm s.t.}, ~ \Vert x \Vert_0\leq s \\
   \text{L0 regularized optimization:} &~\min_{x}&~f(x) +\lambda \Vert x \Vert_0,
   \end{array}
   \end{eqnarray*}
   where $f: \mathbb{R}^{ n}\rightarrow  \mathbb{R}$, $s\ll n, \lambda>0$ and $\Vert x \Vert_0$ is the so-called $\ell_0$ norm that counts the number of nonzero elements of $x$.  <br><br> --->
 
 
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
        <td style="width:9%" align="center"> </td>
        <td style="width:5%" align="center"><a style="font-size: 16px; font-weight: bold; color:#8cd2d5"  href="https://github.com/ShenglongZhou/NHTPver2" target="_blank">NHTP</a></td>
        <td style="width:5%" align="center"><a style="font-size: 16px; font-weight: bold; color:#8cd2d5"  href="https://github.com/ShenglongZhou/NL0R" target="_blank">NL0R</a></td>
        <td style="width:5%" align="center"><a style="font-size: 16px; font-weight: bold; color:#8cd2d5"  href="https://github.com/ShenglongZhou/IIHT" target="_blank">IIHT</a></td>
        <td style="width:5%" align="center"><a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/MIRL1" target="_blank">MIRL1</a></td>
        <td style="width:5%" align="center"><a style="font-size: 16px; font-weight: bold; color:#8cd2d5"  href="https://github.com/ShenglongZhou/HTPCP" target="_blank">HTPCP</a></td>
      </tr>
       <tr>
          <td style="width:9%" align="left"><b>Compressed sensing (CS)</b></td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center"> </td> 
      </tr>
        <tr>
          <td style="width:9%" align="left"><b>Logistic regression (LR)</b></td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center"> </td> 
          <td style="width:5%" align="center"> </td> 
      </tr>
        <tr>
          <td style="width:9%" align="left"><b>Linear complementarity problem (LCP)</b></td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center">$\surd$</td>
          <td style="width:5%" align="center"> </td>
          <td style="width:5%" align="center">$\surd$</td> 
      </tr>
      </table>
 <!---  </details> --->
  </font>



## <b style="font-size:20px">0/1 Loss Optimization Solvers</b>
---

  <font size=4> 
  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/GPSP" target="_blank">GPSP</a> solves the one-bit compressive sensing problems.  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://www.researchgate.net/publication/348371863" target="_blank">Article</a>  
  <p style="line-height: 1;"></p>

  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/Huajun-Wang/L01ADMM" target="_blank">L01ADMM</a> solves the support vector machine. 
  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://arxiv.org/abs/1912.07418" target="_blank">Article</a> 
  </font>
  


## <b style="font-size:20px">EDM Optimization Solvers</b>
---
  
  <font size=4> 

  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/SQREDM" target="_blank">SQREDM</a> solves applications from MDS,  e.g.,  sensor network localization and molecular conformation.  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://ieeexplore.ieee.org/document/8399531" target="_blank">Article</a> 
  <p style="line-height: 1;"></p>

  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://github.com/ShenglongZhou/PREEEDM" target="_blank">PREEEDM</a> solves applications from MDS,  e.g.,  sensor network localization and molecular conformation.   <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://doi.org/10.1007/s12532-019-00168-0" target="_blank">Article</a>  
  </font>




## <b style="font-size:20px">Bilevel Optimization Toolbox </b>
---

  <font size=4>
 
  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://biopt.github.io/bolib/" target="_blank">BOLIBver2</a>, the second version of  the library providing 173 test examples.  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://www.researchgate.net/publication/338375731" target="_blank">Article</a> 
  <p style="line-height: 1;"></p>

  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://biopt.github.io/" target="_blank">BiOpt Toolbox</a>,  bilevel optimization toolbox including <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://biopt.github.io/bolib/" target="_blank">BOLIBver2</a>, 
  <a style="font-size: 16px; font-weight: bold; color:#8cd2d5" href="https://biopt.github.io/solvers/" target="_blank">three solvers</a> and several useful tools. 
  
  </font>
