---
layout: archive
title: ""   
permalink: /SQCQP/
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

 

##  <span style="color:#8C8C8C"> Sparse quadratically constrained quadratic programming solver </span> 
---
Sparse quadratically constrained quadratic programming (SQCQP) takes the form of

\begin{equation}
\begin{array}{rl}
\min_{\mathbf{x}\in\mathbb{R}^{n}} &  \dfrac{1}{2}\mathbf{x}^{\top}\mathbf{Q}_0\mathbf{x}+\mathbf{q}_0^{\top}\mathbf{x}\\\\\\
\mbox{s.t.} & \dfrac{1}{2}\mathbf{x}^{\top}\mathbf{Q}_i\mathbf{x}+\mathbf{q}^{\top}_i\mathbf{x}+c_i\leq0,~i=1,2,\ldots,k\\\\\\
&\mathbf{A}\mathbf{x}\leq \mathbf{b}\\\\\\
&\mathbf{Aeq}\mathbf{x} = \mathbf{beq}\\\\\\
& lb\leq x_i \leq ub,~i=1,2,\ldots,n\\\\\\
&\parallel\mathbf{x}\parallel_0\leq s
\end{array} \tag{SCO}
\end{equation}

where parameters are defined as follows
- $\mathbf{Q}_i\in\mathbb{R}^{n\times n}, \mathbf{q}_i\in\mathbb{R}^{n}, c_i\in\mathbb{R},~~i=0,1,\ldots,k$
- $\mathbf{A}\in\mathbb{R}^{m_1\times n}$, $\mathbf{b}\in\mathbb{R}^{m_1}$
- $\mathbf{Aeq}\in\mathbb{R}^{m_2\times n}$, $\mathbf{beq}\in\mathbb{R}^{m_2}$
- $lb$ and $ub$ are two scalars satisfying $0\in[lb, ub]$
- $\parallel\mathbf{x}\parallel_0$ denotes the so-called $\ell_0$-norm, which counts the number of nonzero entries in $\mathbf{x}$
- $s\ll n$ is a given integer
         
## <span style="color:#8C8C8C"> The solver and its demonstration </span> 
---
<div style="text-align:justify;">
<a style="font-size: 16px; font-weight: bold;color:#006DB0" href="https://github.com/ShenglongZhou/SNSQP" target="_blank">SNSQP</a> was developed based on the algorithm developed in the following paper:
</div>

- <a style="font-size: 14px;color:#000000" href="https://arxiv.org/abs/2503.15109" target="_blank"> S Li, S  Zhou, Z  Luo, Sparse quadratically constrained quadratic programming via semismooth Newton method, *arXiv:2503.15109*, 2025.</a> 

As shown below, inputs need to be specified to call the solver. Note that $\texttt{Qi}$ is a cell that include $Q_i, i=1,2,\ldots,k$ described in (SQCQP).

<p style="line-height: 1;"></p>

```ruby
function Out = SNSQP(n,s,Q0,q0,Qi,qi,ci,ineqA,ineqb,eqA,eqb,lb,ub,pars)

% This code aims at solving the sparse QCQP in the form of
%
%         min             (1/2)(x'{Q_0}x)+q_0'x, 
%         s.t. (1/2)x'*Qi{i}*x+qi(:,i)'*x+ci(i)<=0, i = 1,...,k,
%                                 ineqA*x-ineqb<=0,
%                                      eqA*x-eqb=0,
%                                        lb<=x<=ub,
%                                       ||x||_0<=s,
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
%     Qi:     The quadratic constraint matrix   
%             MUST be a cell array or [], each entry is matrix in R^{n-by-n}          (optional)
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
%             pars.x0       -- Starting point of x                                    (default zeros(n,1))
%             pars.dualquad -- Starting point of mu for quadratic constraint          (default zeros(k,1))
%             pars.dualineq -- Starting point of dual variable for linear inequality  (default zeros(m1,1))
%             pars.dualeq   -- Starting point of dual variable for linear equality    (default zeros(m2,1))
%             pars.dualbd   -- Starting point of nu  for bound/box constraint         (default zeros(n,1))
%             pars.tau      -- A positive scalar                                      (default 1)
%                              NOTE: tuning a proper tau may yield better solutions     
%             pars.itlser   -- Maximum nonumber of line search                        (default 5)
%             pars.itmax    -- Maximum nonumber of iteration                          (default 10000)
%             pars.show     -- Results shown at each iteration if pars.show=1         (default 1)
%                              Results not shown at each iteration if pars.show=0
%             pars.tol      -- Tolerance of the halting condition                    (default 1e-6)
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
Each solver has two required inputs: 'func' defining the example and 'dim' recording dimensions of the example and an optional input 'pars'. Please see the <a style="font-size: 16px; font-weight: bold; color:#007D98" href="\files\menu-of-BiOpt.pdf" target="_blank">menu-of-BiOpt</a> for more details of 'pars'. The chosen solver is <span style="color:#007D98"><b style="font-size:16px">SNLLVF</b></span> and solves one example 'DempeDutta2012Ex24' defined by the following Matlab m-file. This example is from <a style="font-size: 16px; font-weight: bold; color:#007D98"  href="https://biopt.github.io/bolib/" target="_blank">BOLIBver2</a>, where more examples are provided. The <a style="font-size: 16px; font-weight: bold; color:#007D98" href="\files\menu-of-BiOpt.pdf" target="_blank">menu-of-BiOpt</a> also presents other ways to define examples.
</div>

<p style="line-height: 1;"></p>

```ruby
function w=DempeDutta2012Ex24(x,y,keyf,keyxy)
% This file provides all functions defining DempeDutta2012Ex24 problem and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 0 1]
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = (x-1)^2+y^2;
    case 'G'; w = []; 
    case 'f'; w = x^2*y;      
    case 'g'; w = y^2; 
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = 2*(x-1);         
        case 'y' ; w = 2*y;        
        case 'xx'; w = 2;
        case 'xy'; w = 0;
        case 'yy'; w = 2;
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [];    
        case 'y' ; w = [];      
        case 'xx'; w = [];
        case 'xy'; w = [];
        case 'yy'; w = [];
        end           
    case 'f'   
        switch keyxy
        case 'x' ; w = 2*x*y;    
        case 'y' ; w = x^2;          
        case 'xx'; w = 2*y;
        case 'xy'; w = 2*x;
        case 'yy'; w = 0;
        end           
    case 'g'   
        switch keyxy
        case 'x' ; w =   0;  
        case 'y' ; w =   2*y;         
        case 'xx'; w =   0;  
        case 'xy'; w =   0;  
        case 'yy'; w =   2; 
        end        
   end   
end
end
```
 
