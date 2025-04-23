---
layout: archive
title: ""   
permalink: /1BCS/
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

 

##  <span style="color:#8C8C8C"> One-bit compressive sensing</span> 
---
<div style="text-align:justify;">
One-bit compressive sensing (1BCS) problems aim to recover a sparse signal $\mathbf{x}^*\in\mathbb{R}^{n}$ from the following system,
</div>

\begin{equation}
\mathbf{b} = \mathrm{Diag}(\mathbf{h}) \mathrm{sign}(\mathbf{A}\mathbf{x} + \boldsymbol{\varepsilon}) \tag{1bCS}
\end{equation} 

<div style="text-align:justify;">
where $\mathbf{A}\in\mathbb{R}^{m\times n}$ is the sensing matrix, both $\mathbf{b}$ and $\mathbf{h}\in\{-1,1\}^{m}$, $\boldsymbol{\varepsilon}\in\mathbb{R}^{n}$ is the noise, $\mathrm{Diag}(\mathbf{h})$ denotes the diagonal matrix with diagonal entries formed by $\mathbf{h}$, and $\mathrm{sign}(t)$ is the sign function of $t$ defined by $\mathrm{sign}(t)=1$ if $t>0$ and $\mathrm{sign}(t)=-1$ otherwise. Then $\mathrm{sign}(\mathbf{x})=(\mathrm{sign}(x_1),\ldots,\mathrm{sign}(x_n))^\top$. Note that multiplying $\mathrm{Diag}(\mathbf{h})$ means that the sign flips occurs when observed $\mathbf{b}$, making the problem harder. In this model, we assume at most $k$ signs are flipped, namely, $\mathbf{h}$ satisfies $\parallel\mathbf{h}-1\parallel_0\leq k$, where $k$ is a given integer and $\parallel\mathbf{x}\parallel_0$ denotes the $\ell_0$-norm, which counts the number of nonzero entries in $\mathbf{x}$. The following optimization models are solved to recover the signal.
</div> 
 <p style="line-height: 2;"></p>
 <div style="text-align:justify;"> 
◻️ Double-sparsity constrained optimization (DSCO)     
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n},\mathbf{z}\in\mathbb{R}^{m}}~  \parallel  \mathrm{Diag}(\mathbf{b}) \mathbf{A} \mathbf{x}+\mathbf{z} -\epsilon \parallel^2 + \eta \parallel \mathbf{x} \parallel^2,~~~\textrm{s.t.}~ \parallel\mathbf{x} \parallel_0\leq s,~ \parallel \mathbf{z}_+\parallel_0\leq k \tag{DSCO}
\end{equation}
<div style="text-align:justify;">
</div> 
◻️ Step function-regularized optimization (SFRO)  
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n}}~  \sum_{i=1}^n (x_i^2+\varepsilon)^{q/2}+\lambda \parallel (\epsilon- \mathrm{Diag}(\mathbf{b}) \mathbf{A} \mathbf{x})_+ \parallel_0 \tag{SFRO}
\end{equation}
</div> 
<div style="text-align:justify;">
where $s\ll n$, $k\ll m$, $(\epsilon, \eta, \varepsilon, \lambda)>0$, and $\mathbf{z}_+=(\max\{0,z_1\},\ldots,\max\{0,z_m\})^\top$. One can observe that  model (SFRO) is a penalty version of model (DSCO) and term $\parallel \mathbf{z}_+\parallel_0$ is related to the step function (or 0/1 loss function) defined by $\ell_{0/1}(t)=1$ if $t>0$ and $\ell_{0/1}(t)=0$ otherwise, namely, 
  \begin{equation}\|\mathbf{z}_+\|_0=\sum_{i=1}^m \ell_{0/1}\left(z_{i}\right)\nonumber\end{equation}
</div> 
  
---

<div style="text-align:justify;">
The package can be download here - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href="\files\1BCSpack.zip" target="_blank">OBCSpack</a>, which provides 2 solvers from the following papers, where <b style="font-size:16px;color:#777777">GPSP</b> and <b style="font-size:16px;color:#777777">NM01</b> are designed to solve  model (DSCO) and model (SFRO), respectively. 
</div>  

> <b style="font-size:14px;color:#777777">GPSP</b> - <span style="font-size: 14px"> S Zhou, Z Luo, N Xiu, and G Li, Computing one-bit compressive sensing via double-sparsity constrained optimization, IEEE Tran Signal Process, 70:1593-1608, 2022. </span>
<br> <b style="font-size:14px;color:#777777">NM01</b> - <span style="font-size: 14px"> S Zhou, L Pan, N Xiu, and H Qi, Quadratic convergence of smoothing Newton's method for 0/1 loss optimization, SIAM J Optim, 31:3184–3211, 2021. </span>

---

<div style="text-align:justify;">
The inputs and outputs of OBCSpack are detailed below, where inputs $(\texttt{A},\texttt{b},\texttt{s},\texttt{k},\texttt{solver})$ are required. If choose $\texttt{solver=`NM01'}$, then one can set $\texttt{s=[]}$ and $ \texttt{k=[]}$ if they are unkown. The parameters in $\texttt{pars}$ are optional, but setting certain ones may improve the solver's performance and the quality of the solution.
</div>

<p style="line-height: 1;"></p>

```ruby
function out = OBCSpack(A,b,s,k,solver,pars)
% -------------------------------------------------------------------------
% One-bit compressed sensing problem aims to recovere sparse signal x from
%
%                b = Diag(h).*sign( A*x + noise )
%
% 1) The double sparsity constrained optimization (DSCO)
%
%    min  ||Diag(b)*A*x+y-epsilon||^2 + eta||x||^2
%    s.t. ||x||_0<=s, ||y_+||_0<=k
%
% where (epsilon, eta)>0, s\in[1,n], k\in[0,m] are given.
%
% 2) The step function regularized optimization (SFRO)
%
%    min ||x.*x+vareps||^{q/2}_{q/2} + lam*||(epsilon - Diag(b)*A*x)_+||_0
%
% where (vareps, lam, epsilon)> 0, q\in(0,1).  
% -------------------------------------------------------------------------
% Inputs:
%  A:       The sensing matrix \in R^{m-by-n},                   (REQUIRED)
%  b:       The binary observation \in R^m, b_i\in{-1,1}         (REQUIRED)
%  s:       Sparsity level of x, an integer \in[1,n]             (REQUIRED)      
%  k:       An integer in [0,m], e.g., k = ceil(0.01m)           (REQUIRED)       
%  solver:  A text string, can be one of {'GPSP','NM01'}         (REQUIRED)            
%  pars:    Parameters are optional                              (OPTIONAL) 
%           ------------- For GPSP solving (DSCO)--------------------------
%           pars.eps     - The parameter in the model        (default,1e-4)
%           pars.eta     - The penalty parameter       (default,0.01/ln(n))
%           pars.acc     - Acceleration is used if acc=1        (default,0)
%           pars.big     - Start with a bigger s if big=1       (default,1)
%           pars.maxit   - Maximum number of iterations       (default,1e3) 
%           pars.tol     - Tolerance of halting condition    (default,1e-8)
%           -------------  For NM01 solving (SFRO)-------------------------
%           pars.x0      - The initial point           (default,zeros(n,1))
%           pars.q       - Parameter in the objective         (default,0.5)
%           pars.vareps  - Parameter in the objective         (default,0.5)
%           pars.epsilon - Parameter in the objective        (default,0.15)
%           pars.lam     - The penalty parameter                (default,1)
%           pars.tau     - A useful parameter                   (default,1) 
%           pars.maxit   - Maximum number of iterations       (default,1e3)  
% -------------------------------------------------------------------------
% Outputs:
%     out.sol:   The sparse solution x
%     out.time:  CPU time
%     out.iter:  Number of iterations
%     out.obj:   Objective function value at out.sol 
% ------------------------------------------------------------------------
```

<div style="text-align:justify;">
Below is a demonstration of how OBCSpack can be used to solve the problem. You simply need to input the data $(\texttt{A},\texttt{b},\texttt{s},\texttt{k})$ and then choose one solver from $\texttt{\{`GPSP',`NM01'\}}$. 
</div>

<p style="line-height: 1;"></p>

```ruby
clc; close all; clear; addpath(genpath(pwd));

n     = 2000;          % Signal dimension 
m     = ceil(0.5*n);   % Number of measurements
s     = ceil(0.01*n);  % Sparsity level
nf    = 0.05;          % Noisy ratio
r     = 0.02;          % Flipping ratio
k     = ceil(r*m);

A     = randn(m,n);
T     = randperm(n,s);
xo    = zeros(n,1);                      
xo(T) = (1+rand(s,1)).*sign(randn(s,1));  
xo(T) = xo(T)/norm(xo(T));                 % True sparse solution
h     = ones(m,1);                         % Flipping vector
T     = randperm(m,k); 
h(T)  = -h(T);
b     = h.*sign(A(:,T)*xo(T)+nf*randn(m,1));; 

solver = {'GPSP','NM01'};
out    = OBCSpack(A,b,s,k,solver{1});  
fprintf(' Time:                  %6.3f sec\n',out.time);
fprintf(' Absolue error:         %6.2f %%\n', norm(xo-out.sol)*100);
fprintf(' Signal-to-noise ratio: %6.2f\n',-10*log10(norm(xo-out.sol)^2));
fprintf(' Hamming distence:      %6.3f\n',nnz(sign(A*out.sol)-b)/m)
```
