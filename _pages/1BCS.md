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
where $\mathbf{A}\in\mathbb{R}^{m\times n}$ is the sensing matrix, both $\mathbf{b}$ and $\mathbf{h}\in\{-1,1\}^{m}$, $\boldsymbol{\varepsilon}\in\mathbb{R}^{n}$ is the noise, $\mathrm{Diag}(\mathbf{h})$ denotes the diagonal matrix with diagonal entries formed by $\mathbf{h}$, and $\mathrm{sign}(t)$ is the sign function of $t$ defined by $\mathrm{sign}(t)=1$ if $t>0$ and $\mathrm{sign}(t)=-1$ otherwise. Then $\mathrm{sign}(\mathbf{x})=(\mathrm{sign}(x_1),\ldots,\mathrm{sign}(x_n))^\top$. Note that multiplying $\mathrm{Diag}(\mathbf{h})$ means that the sign flips occurs when observed $\mathbf{b}$, making the problem harder. In this model, we assume at most $k$ signs are flipped, namely, $\mathbf{h}$ satisfies $\parallel\mathbf{h}-1\parallel_0\leq k$, where $k$ is a given integer and $\parallel\mathbf{x}\parallel_0$ denotes the $\ell_0$-norm, which counts the number of nonzero entries in $\mathbf{x}$. To recover the signal, the following double-sparsity constrained optimization (DSCO) model is solved, 
</div>      
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{n},\mathbf{y}\in\mathbb{R}^{m}}~  \parallel \mathrm{Diag}(\mathbf{b}) \mathbf{A} \mathbf{x}+\mathbf{y} -\epsilon \parallel^2 + \eta \parallel \mathbf{x} \parallel^2,~~~\textrm{s.t.}~ \parallel\mathbf{x} \parallel_0\leq s,~ \parallel \mathbf{y}_+\parallel_0\leq k \tag{DSCO}
\end{equation}
<div style="text-align:justify;">
where $s\ll n$, $k\ll m$, $\epsilon>0$, $\eta>0$, and $\mathbf{y}_+=(\max\{0,y_1\},\ldots,\max\{0,y_m\})^\top$.
</div> 
  
---
<div style="text-align:justify;">
The solver can be download here - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href="\files\GPSP.zip" target="_blank">GPSP</a>, which was developed from the following paper:
</div>  
>  <span style="font-size: 14px"> S Zhou, Z Luo, N Xiu, and G Li, Computing one-bit compressive sensing via double-sparsity constrained optimization, IEEE Tran Signal Process, 70:1593-1608, 2022. </span>
  
---
<div style="text-align:justify;">
The inputs and outputs of GPSP are detailed below, where $(\texttt{A},\texttt{b},\texttt{s},\texttt{k})$ are required. The parameters in $\texttt{pars}$ are optional, but setting certain ones may improve the solver's performance and the quality of the solution.
</div>

<p style="line-height: 1;"></p>

```ruby
function out = GPSP(A,b,s,k,pars)
% -------------------------------------------------------------------------
% One-bit compressed sensing problem is recovering sparse signal x from
%
%                b = Diag(h).*sign( A*x + noise )
%
% This code aims at solving one-bit compressed sensing 
% via the double sparsity constrained optimization
%
%                min  ||Diag(b)*A*x+y-eps||^2 + eta||x||^2
%                s.t. ||x||_0<=s, ||y_+||_0<=k
%
% where eps>0, eta>0, s\in[1,n], k\in[0,m] are given.
% -------------------------------------------------------------------------
% Inputs:
%     A:    The sensing matrix \in R^{m-by-n},                   (REQUIRED)
%     b:    The binary observation \in R^m, b_i\in{-1,1}         (REQUIRED)
%     s:    Sparsity level of x, an integer \in[1,n]             (REQUIRED)      
%     k:    Upper bound of sign flips of sign(A*x + noise)       (REQUIRED) 
%           An integer in [0,m], e.g., k = ceil(0.01m)         
%     pars: Parameters are all optional                          (OPTIONAL)
%           pars.eps   -- The parameter in the model         (default,1e-4)
%           pars.eta   -- The penalty parameter         (default,.01/ln(n))
%           pars.acc   -- Acceleration is used if acc=1         (default,0)
%           pars.big   -- Start with a bigger s if big=1        (default,1)
%           pars.maxit -- Maximum number of iterations        (default,1e3) 
%           pars.tol   -- Tolerance of halting condition     (default,1e-8)
% Outputs:
%     out.solx : The sparse solution with respect to x in \R^n
%     out.soly : The solution with respect to y in \R^m
%     out.time : CPU time
%     out.iter : Number of iterations
%     out.obj  : Objective function value at (out.solx, out.soly)
% -------------------------------------------------------------------------
```

<div style="text-align:justify;">
Below is a demonstration of how GPSP can be used to solve the problem. You simply need to input the data $(\texttt{A},\texttt{b},\texttt{s},\texttt{k})$. 
</div>

<p style="line-height: 1;"></p>

```ruby
clc; close all; clear;  addpath(genpath(pwd));

n     = 2000;          % Signal dimension 
m     = ceil(0.5*n);   % Number of measurements
s     = ceil(0.01*n);  % Sparsity level
nf    = 0.05;          % Noisy ratio
r     = 0.02;          % Flipping ratio
k     = ceil(r*m);

A     = randn(m,n);
T     = randperm(n,s);
xo    = zeros(n,1);                      
xo(T) = (0.5+rand(s,1)).*sign(randn(s,1));  
xo(T) = xo(T)/norm(xo(T));                 % True sparse solution
bo    = sign(A(:,T)*xo(T)+nf*randn(m,1));
h     = ones(m,1);                         % Flipping vector
T     = randperm(m,k); 
h(T)  = -h(T);
b     = bo.*h; 

out   = GPSP(A,b,s,k); 
fprintf(' Time:                  %6.3f sec\n',out.time);
fprintf(' Absolue error:         %6.3f %%\n', norm(xo-out.solx)*100);
fprintf(' Signal-to-noise ratio: %6.2f\n',-10*log10(norm(xo-out.solx)^2));
fprintf(' Hamming distence:      %6.3f\n',nnz(sign(A*out.solx)-b)/m)
fprintf(' Hamming error:         %6.3f\n',nnz(sign(A*out.solx)-bo)/m)
RecoverShow(xo,out.solx,[1000 450 500 250])
```
