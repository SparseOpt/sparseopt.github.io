---
layout: archive
title: ""   
permalink: /CS/
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

 

##  <span style="color:#8C8C8C"> Compressive sensing</span> 
---
<div style="text-align:justify;">
Compressive sensing (CS) problems aim to recover a signal $\mathbf{x}^*\in\mathbb{R}^{n}$ from the following linear system,
</div>

\begin{equation}
\mathbf{b}=\mathbf{A}\mathbf{x} + \boldsymbol{\epsilon} \tag{CS}
\end{equation}

<div style="text-align:justify;">
where $\mathbf{A}\in\mathbb{R}^{m\times n}$ is the sensing matrix, $\mathbf{b}\in\mathbb{R}^{m}$ is the observation, and $\boldsymbol{\epsilon}\in\mathbb{R}^{n}$ is the noise. 
</div>      

---
<div style="text-align:justify;">
The package can be download here - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href="https://github.com/ShenglongZhou/CSpack" target="_blank">CSpack</a>, which provides 6 solvers from the following papers:</div>  
> <b style="font-size:14px;color:#777777">NHTP</b> - <span style="font-size: 14px"> S Zhou, N Xiu, and H Qi, Global and quadratic convergence of Newton hard-thresholding pursuit, J Mach Learn Res, 22:1-45, 2021. </span>
<br><b style="font-size:14px;color:#777777">GPSP</b> - <span style="font-size: 14px"> S Zhou, Gradient projection Newton pursuit for sparsity constrained optimization, Appl Comput Harmon Anal, 61:75-100, 2022. </span>
<br><b style="font-size:14px;color:#777777">PSNP</b> - <span style="font-size: 14px"> S Zhou, X Xiu, Y Wang, and D Peng, Revisiting Lq ( 0 ≤ q < 1 ) norm regularized optimization, arXiv:2306.14394, 2023. </span>
<br><b style="font-size:14px;color:#777777">NL0R</b> - <span style="font-size: 14px"> S Zhou, L Pan, and N Xiu, Newton method for l0 regularized optimization, Numer Algorithms, 88:1541–1570, 2021. </span>
<br><b style="font-size:14px;color:#777777">IIHT</b> - <span style="font-size: 14px"> L Pan, S Zhou, N Xiu, and H Qi, A convergent iterative hard thresholding for nonnegative sparsity optimization, Pac J Optim, 13:325-353, 2017. </span>
<br><b style="font-size:14px;color:#777777">MIRL1</b> - <span style="font-size: 14px"> S Zhou, N Xiu, et. al, A Null-space-based weighted l1 minimization approach to compressed sensing, Inf inference, 5:76-102, 2016. </span>

---
<div style="text-align:justify;">
As shown below, the solver inputs are described in details. Note that input $\texttt{Qi}$ is a cell that include $Q_i, i=1,2,\ldots,k$ described in (SQCQP). If some constraints are absent, then just put them as an empty set, i.e.,  $\texttt{[ ]}$.
</div>

<p style="line-height: 1;"></p>

```ruby
function Out = SNSQP(n,s,Q0,q0,Qi,qi,ci,ineqA,ineqb,eqA,eqb,lb,ub,pars)

% This code aims at solving the sparse QCQP in the form of
%
%         min  (1/2)(x'{Q_0}x)+q_0'x  
%         s.t. (1/2)x'*Qi{i}*x+qi(:,i)'*x+ci(i)<=0, i = 1,...,k 
%                                 ineqA*x-ineqb<=0 
%                                      eqA*x-eqb=0 
%                                        lb<=x<=ub 
%                                       ||x||_0<=s 
% where Qi = {Qi{1},...,Qi{k}}, Qi{i} \in R^{n-by-n}, qi \in R^{n-by-k},  ci \in R^{k}
%       ineqA \in R^{m1-by-n},  ineqb \in R^{m1} 
%       eqA   \in R^{m2-by-n},  eqb   \in R^{m2}
%       s << n

%---------------------------------------------------------------------------------------------------           
% Inputs:
%     n:      Dimension of the solution x                                             (required)
%     s:      Sparsity level of x, an integer between 1 and n-1                       (required)
%     Q0:     The quadratic objective matrix in R^{n-by-n}                            (required)        
%     q0:     The quadratic objective vector in R^n                                   (required)
%     Qi:     The quadratic constraint matrix                                         (optional) 
%             MUST be a cell array or [], each entry is matrix in R^{n-by-n}           
%     qi:     The quadratic constraint vector. MUST be a matrix in R^{n-by-k} or []   (optional)           
%     ci:     The quadratic constraint constant in R, must be a vector or []          (optional)
%     ineqA:  The linear inequality constraint matrix in R^{m1-by-n}   or []          (optional)
%     ineqb:  The linear inequality constraint vector in R^{m1}        or []          (optional)
%     eqA:    The linear equality constraint matrix in R^{m2-by-n}     or []          (optional)
%     eqb:    The linear equality constraint vector in R^{m2}          or []          (optional)
%     lb:     The lower bound of x                                                    (optional)
%     ub:     The upper bound of x                                                    (optional)
%             NOTE: 0 must in [lb ub]
%     pars:   Parameters are all OPTIONAL
%             pars.x0       -- Initial point of x                                     (default zeros(n,1))
%             pars.dualquad -- Initial point of mu for quadratic constraints          (default zeros(k,1))
%             pars.dualineq -- Initial point of dual variable for linear inequalities (default zeros(m1,1))
%             pars.dualeq   -- Initial point of dual variable for linear equalities   (default zeros(m2,1))
%             pars.dualbd   -- Initial point of nu  for bound/box constraints         (default zeros(n,1))
%             pars.tau      -- A positive scalar                                      (default 1)
%                              NOTE: tuning a proper tau may yield better solutions     
%             pars.itlser   -- Maximum nonumber of line search                        (default 5)
%             pars.itmax    -- Maximum nonumber of iteration                          (default 10000)
%             pars.show     -- Results shown at each iteration if pars.show=1         (default 1)
%                              Results not shown at each iteration if pars.show=0
%             pars.tol      -- Tolerance of the halting condition                     (default 1e-6)
%
% Outputs:
%     Out.sol:           The sparse solution x
%     Out.sparsity:      Sparsity level of Out.sol
%     Out.error:         Error used to terminate this solver
%     Out.time           CPU time
%     Out.iter:          Number of iterations
%     Out.obj:           Objective function value at Out.sol
%---------------------------------------------------------------------------------------------------
```

<div style="text-align:justify;">
Below is a demonstration of the solver applied to a sparse portfolio selection problem. The parameters in $\texttt{pars}$ are optional; however, specifying certain ones can enhance the solver's performance and solution quality.
</div>

<p style="line-height: 1;"></p>

```ruby
% demon sparse portfolio selection problems
clc; clear all; close all;  addpath(genpath(pwd));

n     = 1000;
s     = 10;

B     = 0.01 * rand(ceil(n/4),n);
D     = diag(0.01*rand(n,1));
Q0    = 2*( B'*B + D);
q0    = zeros(n,1); 
Qi    = cell(1,1);
Qi{1} = 2*D;
qi    = zeros(n,1);
ci    = -0.001;
ineqA = -0.5*randn(1,n);
ineqb = -0.002;
eqA   = ones(1,n);
eqb   = 1;
lb    = 0;
ub    = 0.3;
    
pars.x0       = ((lb+ub)/2).*ones(n,1);
pars.tau      = 1; % decrease this value if the algorithm do not converge
pars.dualquad = 0*ones(length(ci));
pars.dualineq = 0.001*ones(length(ineqb)); 
pars.dualeq   = 0.001*ones(length(eqb));
Out           = SNSQP(n,s,Q0,q0,Qi,qi,ci,ineqA,ineqb,eqA,eqb,lb,ub,pars);
```
 
