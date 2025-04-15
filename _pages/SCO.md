---
layout: archive
title: ""   
permalink: /solvers/
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

## <span style="color:#225675"><b style="font-size:20px">Bilevel optimization solvers</b></span> 
---
Bilevel optimization has the form

\begin{eqnarray}\min_{x,y} &&   F(x,y) \nonumber\\\\\\
\mbox{s.t.} && G(x,y)\leq 0,~ y\in \mbox{argmin}_y~ \lbrace f(x,y)\mid g(x,y)\leq 0\rbrace, \nonumber
\end{eqnarray}
<div style="text-align:justify;">
where  $F,f:\mathbb{R}^{n_x}\times\mathbb{R}^{n_y}\rightarrow \mathbb{R}$, $G:\mathbb{R}^{n_x}\times\mathbb{R}^{n_y}\rightarrow \mathbb{R}^{n_G}$ and $g:\mathbb{R}^{n_x}\times\mathbb{R}^{n_y}\rightarrow \mathbb{R}^{n_g}$.  A problem is linear if all involved functions are linear and nonlinear otherwise. We call a problem <i>a simple bilevel optimization</i> if there is no variable $x$ involved. This two level optimization can be transformed into a single-level version so that Semi-smooth Newton type method is able to be employed. 
</div>
 

### <span style="color:#225675"><b style="font-size:20px">Bilevel optimization reformulations</b></span> 
---

<div style="text-align:justify;">
<span style="color:#007D98"><b style="font-size:16px">SNLLVF:</b></span> This approach is to transform the bilevel program into a single-level optimization problem by using
the lower-level value function (LLVF) reformulation, namely,  $g(x,y)\leq 0, f(x,y) = \varphi(x) := \underset{z}\min\lbrace f(x,z) \mid g(x,z)\leq 0 \rbrace$. By doing so,  SNLLVF  aims at solving a partial penalization
</div>
\begin{eqnarray}\min_{x,y} && F(x,y) + \lambda (f(x,y) -\varphi(x)) \nonumber\\\\\\
\mbox{ s.t. } &&  G(x,y)\leq 0,~ g(x,y)\leq 0. \nonumber
\end{eqnarray}
 
 
<div style="text-align:justify;">
<span style="color:#007D98"><b style="font-size:16px">SNQVI:</b></span> This approach is to transform the bilevel program into a single-level optimization problem 
by converting the lower-level problem to the quasi-variational inequality conditions: $ \langle \nabla_y f(x,y), z-y \rangle \geq0, \forall~z: g(x,z)\leq 0$ for any $y: g(x,y)\leq 0$. Let $\varphi(x,y) := \underset{z}\min\lbrace\nabla_y f(x,y)^\top z \mid g(x,z)\leq 0  \rbrace$.  By doing so,  SNQVI  aims at solving a partial penalization
</div>
\begin{eqnarray}\min_{x,y} && F(x,y)+ \lambda ( \nabla_y f(x,y)^\top y-\varphi(x,y) ) \nonumber\\\\\\
\mbox{ s.t. } && G(x,y)\leq 0,  \ \   g(x,y)\leq 0. \nonumber
\end{eqnarray}


 
<div style="text-align:justify;">
<span style="color:#007D98"><b style="font-size:16px">SNKKT:</b></span> This approach is to transform the bilevel program into a single-level optimization problem 
by converting the lower-level problem to the KKT conditions: $ \nabla_y f(x,y)-\nabla_y g(x,y)^\top z=0,~ g(x,y)\leq 0,~  z \leq 0,~   g(x,y)^\top z=0. $ By doing so, SNKKT   aims at solving a partial penalization
 </div>
\begin{eqnarray}
\min_{x,y,z}  && F(x,y)+ \lambda g(x,y)^\top z\nonumber\\\\\\
\mbox{ s.t. } && G(x,y)\leq 0,  \   g(x,y)\leq 0,\ z \leq 0, \ \nabla_y f(x,y)-\nabla_y g(x,y)^\top z=0. \nonumber
\end{eqnarray}  


### <span style="color:#225675"><b style="font-size:20px">Demonstration of solvers</b></span> 
---
<div style="text-align:justify;">
<a style="font-size: 16px; font-weight: bold; color:#007D98" href="https://biopt.github.io/solvers/" target="_blank">BiOpt-Solvers</a> provides three solvers:  <span style="color:#007D98"><b style="font-size:16px">SNLLVF</b></span>, <span style="color:#007D98"><b style="font-size:16px">SNQVI</b></span> and  <span style="color:#007D98"><b style="font-size:16px">SNKKT</b></span> based on the above three reformulations. Detailed instructions can be found in  the <a style="font-size: 16px; font-weight: bold; color:#007D98" href="\files\menu-of-BiOpt.pdf" target="_blank">menu-of-BiOpt</a>. Here we give a simple example to illustrate it:
</div> 

<p style="line-height: 1;"></p>

```
clc; clear; close all; 

ExName     = 'DempeDutta2012Ex24'; 
func       = str2func(ExName);
dim        = [1 1 0 1];     % required
pars.xy    = [1;1];         % optional

pars.lam   = 1;             % optional
pars.keep  = 0;             % optional 
pars.check = 1;             % optional

SolNo      = 1;     % choose a solver
Solvers    = {'SNLLVF','SNQVI','SNKKT'}; 
solver     = str2func(Solvers{SolNo});  
Out1       = solver(func, dim,  pars);
```

<div style="text-align:justify;">
Each solver has two required inputs: 'func' defining the example and 'dim' recording dimensions of the example and an optional input 'pars'. Please see the <a style="font-size: 16px; font-weight: bold; color:#007D98" href="\files\menu-of-BiOpt.pdf" target="_blank">menu-of-BiOpt</a> for more details of 'pars'. The chosen solver is <span style="color:#007D98"><b style="font-size:16px">SNLLVF</b></span> and solves one example 'DempeDutta2012Ex24' defined by the following Matlab m-file. This example is from <a style="font-size: 16px; font-weight: bold; color:#007D98"  href="https://biopt.github.io/bolib/" target="_blank">BOLIBver2</a>, where more examples are provided. The <a style="font-size: 16px; font-weight: bold; color:#007D98" href="\files\menu-of-BiOpt.pdf" target="_blank">menu-of-BiOpt</a> also presents other ways to define examples.
</div>

<p style="line-height: 1;"></p>

```
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
 
