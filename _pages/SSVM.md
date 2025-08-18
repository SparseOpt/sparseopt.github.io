---
layout: archive
title: ""   
permalink: /SSVM/
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

 

##  <span style="color:#8C8C8C"> Sparse support vector machine</span> 
---
<div style="text-align:justify;">
  The  soft-margin support vector machine (SVM) takes the form of 
</div>

\begin{equation}
\min_{(\mathbf{w};b)\in\mathbb{R}^{n+1}}~\frac{1}{2}\parallel \mathbf{w} \parallel^2 + C \sum_{i=1}^m\ell\left(1-y_i(b+ \mathbf{a}_i^\top\mathbf{w})\right) \tag{SVM}
\end{equation} 

<div style="text-align:justify;">
where $\mathbf{A}=(\mathbf{a}_1,\ldots,\mathbf{a}_m)^\top\in\mathbb{R}^{m\times n}$ is the sample matrix, $\mathbf{y}=(y_1,\ldots,y_m)^\top\in\mathbb{R}^{m}$ is the label vector, $C>0$ is the penalty parameter, and $\ell$ is a loss function. One popular loss function is the hinge loss defined by  $\ell_{h}(t)=\max\{t,0\}.$ Two loss functions are considered below, resulting in two SVM models.
</div>      

<p style="line-height: 2;"></p>

◻️ Step function regularized SVM
\begin{equation}
\min_{(\mathbf{w};b)\in\mathbb{R}^{n+1}}~\frac{1}{2}\parallel \mathbf{w} \parallel^2 + C \sum_{i=1}^m\mathrm{step}\left(1-y_i(b+  \mathbf{a}_i^\top\mathbf{w})\right) \tag{L01SVM}
\end{equation} 
<div style="text-align:justify;">
where $\mathrm{step}(t)$ is the step function (or 0/1 loss function) defined by $\mathrm{step}(t)=1$ if $t>0$ and $\mathrm{step}(t)=0$ otherwise. By letting $\parallel\mathbf{x}\parallel_0$ denote the L0 norm (i.e., the number of nonzero entries) of $\mathbf{x}$ and $\mathbf{z}_+=(\max\{0,z_1\},\ldots,\max\{0,z_m\})^\top$, it follows $\sum_{i=1}^m\mathrm{step}\left(1-y_i(b+  \mathbf{a}_i^\top\mathbf{w})\right)$=$\| (1-\mathbf{A}\mathbf{w}-b\mathbf{y} )_+ \|_0$. 
</div>

<!--
◻️ $\ell_{cC}$ regularized  SVM
\begin{equation}
\min_{(\mathbf{w};b)\in\mathbb{R}^{n+1}}~\frac{1}{2} \parallel  \mathbf{w} \parallel^2 + \sum_{i=1}^m\ell_{cC}\left(1-y_i(b+  \mathbf{a}_i^\top\mathbf{w})\right) \tag{SFRSVM}
\end{equation} 
<div style="text-align:justify;">
where  $\ell_{cC}(t)=Ct^2/2$ if $t>0$ and $\ell_{cC}(t)=ct^2/2$ otherwise with $C>c>0$. The dual problem of (LcCSVM) is the following quadratic kernel-based SVM problem
</div>  

\begin{equation}
\min_{\boldsymbol{\alpha}\in\mathbb{R}^{m}}~\frac{1}{2} \boldsymbol{\alpha}^\top \mathbf{Q} \boldsymbol{\alpha} +\frac{1}{2}\sum_{i=1}^m h_{cC}(\alpha_i) -\mathbf{e}^\top\boldsymbol{\alpha}, ~~~~ \text{s.t.} ~~\mathbf{y}^\top\boldsymbol{\alpha}=0\tag{QKSVM}
\end{equation} 
<div style="text-align:justify;">
where $\mathbf{Q}=(Q_{ij})_{1\leq i,j\leq m}$ with $Q_{ij}=y_iy_j\mathbf{a}_i^\top\mathbf{a}_j$, $\mathbf{e}=(1,\ldots,1)^\top$, and $h_{cC}(t)=t^2/C$ if $t>0$ and $\ell_{cC}(t)=t^2/c$.
</div>  
-->

◻️ Sparsity constrained quadratic kernel-based SVM 
\begin{equation}
\min_{\boldsymbol{\alpha}\in\mathbb{R}^{m}}~\frac{1}{2} \boldsymbol{\alpha}^\top \mathbf{Q} \boldsymbol{\alpha} +\frac{1}{2}\sum_{i=1}^m h_{cC}(\alpha_i) -\mathbf{e}^\top\boldsymbol{\alpha}, ~~~~ \text{s.t.} ~~\mathbf{y}^\top\boldsymbol{\alpha}=0,~\parallel  \boldsymbol{\alpha} \parallel_0\leq s \tag{SCSVM}
\end{equation} 
<div style="text-align:justify;">
where $\mathbf{Q}\in\mathbb{R}^{m\times m}$ with each entry $Q_{ij}=y_iy_j\mathbf{a}_i^\top\mathbf{a}_j$, $\mathbf{e}=(1,\ldots,1)^\top$, and $h_{cC}(t)=t^2/C$ if $t>0$ and $\ell_{cC}(t)=t^2/c$,  $C>c>0$, and $s\ll m$. In fact, model (SCSVM) without sparsity constraint $\parallel  \boldsymbol{\alpha} \parallel_0\leq s$ is the dual problem of model (SVM) with $\ell=\ell_{cC}$, where  $\ell_{cC}(t)=t^2/2$ if $t>0$ and $\ell_{cC}(t)=(c/C)t^2/2$ otherwise. 
</div>  

> <div style="text-align:justify;"> According to the Representer Theorem,  optimal solution $ \mathbf{w}^* $ to (SVM) and optimal solution $\boldsymbol{\alpha}^* $ to the dual SVM satisfy $ \mathbf{w}^* = \sum_{i=1}^m \alpha_i^* y_i \mathbf{a}_i $. The training vectors $ \mathbf{a}_i $ corresponding to nonzero $ \alpha_i^* $ are known as support vectors. Therefore, both models (SFRSVM) and (SCSVM) enable the reduction of support vectors. </div> 

---
<div style="text-align:justify;">
The package can be downloaded here - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href="\files\SSVMpack.zip" target="_blank">SSVMpack</a>, which provides 2 solvers from the following 2 papers, where <b style="font-size:16px;color:#777777">NM01</b> and <b style="font-size:16px;color:#777777">NSSVM</b> are designed to solve (SFRSVM) and (SCSNM), respectively.
</div>  

> <div style="text-align:justify;"> <b style="font-size:14px;color:#777777">NM01</b> -<span style="font-size: 14px"> S Zhou, L Pan, N Xiu, and H Qi, Quadratic convergence of smoothing Newton's method for 0/1 loss optimization, SIOPT, 31:3184-3211, 2021. </span> </div>
> <div style="text-align:justify;">  <b style="font-size:14px;color:#777777">NSSVM</b> -<span style="font-size: 14px"> S Zhou, Sparse SVM for sufficient data reduction, IEEE TPAMI, 44:5560-5571, 2022. </span> </div>

---
<div style="text-align:justify;">
Below is a demonstration of how <b style="font-size:16px;color:#777777">SSVMpack</b> can be used to solve the problem. You simply need to input the data ($\texttt{A}$, $\texttt{y}$) and select $\texttt{solver}$ from {'$\texttt{NM01}$', '$\texttt{NSSVM}$'}. The parameters in $\texttt{pars}$ are optional, but setting certain ones (particularly, $\texttt{pars.C}$ and $\texttt{pars.s0}$)  can improve the solver's performance and the quality of the solution.
</div>

<p style="line-height: 1;"></p>

```ruby
clc; close all; clear all; addpath(genpath(pwd));
  
load dhrb.mat; 
load dhrbclass.mat;  
[m0,n]  = size(A);         
A       = normalization(A,2*(max(A(:))>1)); % normalize the data
% Split the data into training and testing sets
m       = ceil(0.9*m0);         
Train   = randperm(m0,m);  
Ttest   = setdiff(1:m0,Train);  
Atrain  = A(Train,:);  
ytrain  = y(Train,:); 
Atest   = A(Ttest,:);  
ytest   = y(Ttest);    

pars.C  = 0.25;
solver  = {'NM01','NSSVM'};
out     = SSVMpack(Atrain,ytrain,solver{1},pars);
acc     = accuracy(Atrain,out.w,ytrain);
tacc    = accuracy(Atest,out.w,ytest);

fprintf(' Training  Time:             %5.3fsec\n',out.time);
fprintf(' Training  Size:             %dx%d\n',size(Atrain,1),size(Atrain,2));
fprintf(' Training  Accuracy:         %5.2f%%\n', acc*100);
fprintf(' Testing   Size:             %dx%d\n',size(Atest,1),size(Atest,2));
fprintf(' Testing   Accuracy:         %5.2f%%\n',tacc*100);
fprintf(' Number of Support Vectors:  %d\n',out.sv); 
```
<div style="text-align:justify;">
The citation for <b style="font-size:16px;color:#777777">SSVMpack</b> is shown below. Inputs ($\texttt{A}$, $\texttt{y}$, $\texttt{solver}$) are required, $\texttt{solver}$ is chosen from {'$\texttt{NM01}$', '$\texttt{NSSVM}$'}.   If $\texttt{solver}$='$\texttt{NSSVM}$', then set a proper $\texttt{pars.s0}$ can enhance solution quality.  Another important parameter is $\texttt{pars.C}$, which can be tuned using the Cross-validation.
</div>

<p style="line-height: 1;"></p>

```ruby
function out = SSVMpack(A,y,solver,pars)
% -------------------------------------------------------------------------
% This package aims to solve the binary classification problems
% Inputs:
%  A:       The smaple matrix \in R^{m-by-n},                    (REQUIRED)
%  y:       The binary label \in R^m, b_i\in{-1,1}               (REQUIRED)    
%  solver:  A text string, can be one of {'NM01','NSSVM'}        (REQUIRED)            
%  pars:    Parameters are optional                              (OPTIONAL) 
%           -------------  For NSSVM --------------------------------------
%           pars.alpha --  Starting point in \R^m       (default zeros(m,1))
%           pars.s0    --  The initial sparsity    (default n(log(m/n))^2))
%           pars.C     --  A positive scalar in (0,1]         (default 1/4)  
%           pars.c     --  A positive scalar in (0,1]         (default 1/8)  
%           pars.tune  --  Tune the sparsity level              
%                          Do not tune the sparsity level       (default 0)
%           pars.maxit --  Maximum number of iterations      (default e000) 
%           pars.tol   --  Tolerance of the halting criteria (default 1e-4) 
%           pars.disp  --  Display results for each step        (default 1)  
%                          Do not display results for each step 
%           -------------  NM01 -------------------------------------------
%           pars.x0    --  The initial point           (default zeros(n,1))
%           pars.C     --  The penalty parameter                (default 1)
%           pars.tau   --  A useful paramter                    (default 5)
%           pars.maxit --  Maximum number of iterations      (default 1000)  
%           pars.tol   --  Tolerance of the halting criteria (default 1e-4) 
%           pars.disp  --  Display results for each step        (default 1)  
%                          Do not display results for each step 
% -------------------------------------------------------------------------
% Outputs:
%     out.w:      The solution of the primal problem, i.e., the classifier
%     out.sv:     Number of support vectors 
%     out.time    CPU time
%     out.iter:   Number of iterations
%     Out.acc:    Classification accuracy
% -------------------------------------------------------------------------
% Send your comments and suggestions to <slzhou2021@163.com> 
% Warning: Accuracy may not be guaranteed !!!!!! 
% -------------------------------------------------------------------------
```
