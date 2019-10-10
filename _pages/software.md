---
layout: archive
title: ""  
permalink: /software/
author_profile: true
redirect_from:
  - /resume
---

<span style="color:grey">Sparse optimization solvers</span> 
---

<span style="color:orange"> {\bf Solvers}</span> 

* [NHTP](https://github.com/ShenglongZhou/NHTP), a matlab package aiming at solving the sparsity constrained optimization problems including compressed sensing, sparse logistic regression.  This is the code source for paper [Global and Quadratic Convergence of Newton Hard-Thresholding Pursuit](https://shenglongzhou.github.io/publication/2019-01-01-Global-and-Quadratic-Convergence-of-Newton-Hard-Thresholding-Pursuit).

* [IIHT](https://github.com/ShenglongZhou/IIHT), a matlab package aiming at solving the sparsity constrained optimization problems . This is the code source for paper [ A Convergent Iterative Hard Thresholding for Sparsity and Nonnegativity Constrained Optimization](http://www.ybook.co.jp/online2/oppjo/vol13/p325.html). 

* [MIRL1](https://github.com/ShenglongZhou/MIRL1), a matlab package aiming at solving the reweighted $\ell_1$ minimization. This is the code source for paper [A Null-space-based Weighted $\ell_1$ Minimisation Approach to Compressed Sensing](https://shenglongzhou.github.io/publication/2016-01-01-A-Null-space-based-Weighted-l1-Minimisation-Approach-to-Compressed-Sensing).
 
 
* [HTPCP](https://github.com/ShenglongZhou/HTPCP), a matlab package aiming at solving the sparse linear/nonlinear complementarity problems . This is the code source for paper [A Half Thresholding Projection Algorithmfor Sparse Solutions of LCPs](https://link.springer.com/article/10.1007/s11590-014-0834-7). 
 
<span style="color:orange"> {\bf  Applications</span> 
 
 * Compressed sensing:
 \begin{eqnarray}CSS: && \min_{x} ~ 0.5\Vert Ax-b \Vert^2, ~ {\rm s.t.}, ~ \Vert x\Vert _0\leq s  \nonumber\\\\\\
CSL_1:&& \min_{x} ~ 0.5\Vert Wx\Vert_1, ~ \mbox{s.t.},~ Ax=b. \nonumber
\end{eqnarray}

* logistic regression:
\begin{eqnarray}
\min_{x}~  \frac{1}{m}\sum_{i=1}^{m}\left\lbrace \ln(1+ e^{\langle a_i, x\rangle})-b_i\langle a_i, x\rangle\right \rbrace+\mu\Vert x\Vert_2^2 , {\rm s.t.},\ \Vert x\Vert_0\leq s.
\end{eqnarray}


<span style="color:grey">Bilevel optimization toolbox</span> 
---

* [BOLIB](https://github.com/ShenglongZhou/BOLIB), the first version of the library providing 124 test examples. See [BOLIB: bilevel optimization library of test problems](https://arxiv.org/abs/1812.00230) for more information.

* [BOLIBver2](https://biopt.github.io/bolib/), the second version of  the library providing 173 test examples. See [BOLIB2019: bilevel optimization library of test problems version 2](https://biopt.github.io/bolib/) for more information.

* [BiOpt Toolbox](https://biopt.github.io/),  bilevel optimization toolbox including 173 test examples from [BOLIBver2](https://biopt.github.io/bolib/), three solvers to solve bilevel optimization problems and several useful tools.

<span style="color:grey">Euclidean distance matrix optimization solvers</span> 
---

* [SQREDM](https://github.com/ShenglongZhou/SQREDM), a matlab package aiming at solving Euclidean distance matrix optimization problems from multidimensional scaling including sensor network localization, molecular conformation.  This is the code source for paper [A Fast Matrix Majorization-Projection Method for Penalized Stress Minimization with Box Constraints](https://shenglongzhou.github.io/publication/2018-07-28-A-Fast-Matrix-Majorization-Projection-Method-for-Penalized-Stress-Minimization-with-Box-Constraints).

* [PREEEDM](https://github.com/ShenglongZhou/SQREDM), a matlab package aiming at solving Euclidean distance matrix optimization problems from multidimensional scaling including sensor network localization, molecular conformation.  This is the code source for paper [Robust Euclidean Embedding via EDM Optimization](https://shenglongzhou.github.io/publication/2019-08-09-Robust-Euclidean-Embedding-via-EDM-Optimisation).

