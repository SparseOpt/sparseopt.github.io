---
layout: archive
title: ""   
permalink: /SFCO/
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

## Step function-constrained optimization
---

<p style="line-height: 1;"></p>
\begin{equation}
\min_{\mathbf{x}\in\mathbb{R}^{K}} ~~  f(\mathbf{x}),~~~~ \mbox{s.t.}~~ \parallel\mathbf{G}(\mathbf{x})\parallel_0^+\leq s,~\mathbf{x}\in \Omega  \tag{SFCO}
\end{equation}

<div style="text-align:justify;">
where $\mathbf{G}(\mathbf{x})\in\mathbb{R}^{M \times N}$ with each the ($i,j$)th entry being $G_{ij}(\mathbf{x})$, functions $f:\mathbb{R}^{K}\rightarrow \mathbb{R}$ and $G_{ij}:\mathbb{R}^{K}\rightarrow \mathbb{R}$ are (preferably twice) continuously differentiable, $s\ll N$ is an integer, and $\Omega\subseteq\mathbb{R}^{K}$ is a closed and convex set. For a matrix $\mathbf{Z}\in\mathbb{R}^{M \times N}$,  measure $\|\mathbf{Z}\|_0^+$ counts the number of its columns that have positive entries, namely, 
  \begin{equation}\|\mathbf{Z}\|_0^+= \mathrm{step}\Big(\max_{i=1,\ldots,M} Z_{i1}\Big)+\cdots+\mathrm{step}\Big(\max_{i=1,\ldots,M} Z_{iN}\Big)\nonumber\end{equation}
  Here, $\mathrm{step}(t)$ is the step function (or 0/1 loss function) defined by $\mathrm{step}(t)=1$ if $t>0$ and $\mathrm{step}(t)=0$ otherwise. When $M=1$, matrix $\mathbf{Z}$ reduces to a vector $\mathbf{z}\in\mathbb{R}^{N}$, let $\mathbf{z}_+=$ $(\max\{0,z_1\}$ $\ldots$ $\max\{0,z_N\})^\top$ and  $\parallel\mathbf{z}\parallel_0$ denote its L0 norm that counts the number of its nonzero entries. As a result,  
  \begin{equation*}\|\mathbf{z}\|_0^+= \mathrm{step}(z_1)+\cdots+\mathrm{step}(z_N)=\|\mathbf{z}_+\|_0\end{equation*}
In the present package, $\Omega$ can be one of the following sets:
</div>  

- A sphare: $\lbrace\mathbf{x}: \parallel\mathbf{x}\parallel\leq r\rbrace$, where $r>0$: 
- A halfspace: $\lbrace\mathbf{x}: \mathbf{a}^T\mathbf{x}\leq b\rbrace$, where $\mathbf{a}\in\mathbb{R}^{K}$ and $b\in\mathbb{R}$
- A hyperplane: $\lbrace\mathbf{x}: \mathbf{A} \mathbf{x}=  \mathbf{b}\rbrace$, where $\mathbf{A}\in\mathbb{R}^{S\times K}$ and $ \mathbf{b}\in\mathbb{R}^{S}$
- A box:  $\lbrace\mathbf{x}: l\leq x_i \leq u, i=1,\ldots,K\rbrace$, where  $l \leq u$ can be $-\infty$ and $+\infty$. Hence, when $l=-\infty$ and $u=+\infty$, set $\Omega=\mathbb{R}^{K}$, which means no such constraint;  when $l=0$ and $u=+\infty$, set $\Omega$ is a non-negative cone


---
<div style="text-align:justify;"> 
Package - <a style="font-size: 16px; font-weight: bold;color:#006DB0" href=" " target="_blank">SFCOpack-Matlab</a> or <a style="font-size: 16px; font-weight: bold;color:#006DB0" href=" " target="_blank">SFCOpack-Python</a>  provides 1 solver from the following paper:
</div>

> <b style="font-size:14px;color:#777777">SNSCO</b> - <span style="font-size: 14px"> S Zhou, L Pan, N Xiu,  and G  Li, A 0/1 constrained optimization solving sample average approximation for chance constrained programming, MOR, 2024. </span>

---
<div style="text-align:justify;">  
Note that $\texttt{SNSCO}$ is a second-order method, requiring the objective, gradient, and Hessian of functions $f(\mathbf{x})$ and $\mathbf{G}(\mathbf{x})$. Based on Matlab syntax (similar to Python syntax), below is an example of how to define these functions for the solver to solve a recovery problem. The problem has $f(\mathbf{x})=0.5\parallel\mathbf{B}\mathbf{x}-\mathbf{d}\parallel^2$ and $\mathbf{G}_{ij}(\mathbf{x})= \langle\mathbf{A}_{:(j-1)M+i}, \mathbf{x}\rangle-C_{ij}$. The two MATLAB code snippets below respectively define the function values, gradients, and Hessians of $f$ and $\mathbf{G}$. For example, in the inputs of the function handle $\texttt{FuncfRecovery}$, $\texttt{x}$ is the variable, while（$\texttt{B}$, $\texttt{d}$, $\texttt{BtB}$）are the data involved in the function $f(\mathbf{x})$. When calling the function $\texttt{FuncfRecovery}$, users need to define these data.
</div>
<p style="line-height: 1;"></p>

```ruby
function [objef, gradf, hessf] = FuncfRecovery(x,B,d,BtB)
    % This code provides information for an objective function
    %     f(x) = 0.5*||Bx-d||^2  
    % x is the variable 
    % (B,d,BtB) are data and need to be input
    % B\in R^{m by n}, d\in R^{m by 1}, and BtB = B'*B
    
    Bxd   = B*x-d;
    objef = norm(Bxd)^2/2;  % objective
    gradf = (Bxd'*B)';      % gradient
    hessf = BtB;            % Hessian
    clear B d BtB
end
```

```ruby
function [G,gradG,gradGW, hessGW] = FuncGRecovery(x,W,Ind,A,C,K,M,N) 
    % This code provides information for function G(x): R^K -> R^{M x N}
    % (x,W,Ind) are variables
    % (A,C,K,M,N) are data and parameters, and need to be input 
    % For each i=1,...,M and j=1,...,N:
    %           G_{ij}(x) = <A_{:,(j-1)*M+i}, x>^2 - C_{ij} 
    % where A_{:,(j-1)*M+i} is the ((j-1)*M+i)th column of A
    % A \in R^{K by M*N} and C \in R^{M by N} 
    
    G0  = x'*A;
    G   = reshape(G0,M,N);
    G   = G.^2-C;                                         % function
    if  nargout > 1
        if  isempty(Ind) 
            gradG   = [];
            gradGW  = zeros(K,1);                         % gradient
            hessGW  = zeros(K,1);                         % Hessian
        else 
            AInd    = A(:,Ind);    
            gradG   = AInd.*G0(Ind); 
            WInd    = W(Ind);
            gradGW  = gradG*reshape(WInd,length(Ind),1);  % gradient
            hessGW  = AInd*diag(WInd)*AInd.';             % Hessian     
        end
    end  
    clear A C K M N
end
```

<div style="text-align:justify;">  
After defining functions $f$ and $\mathbf{G}$, users need to select constraint set $\Omega$. Currently, $\Omega$ can be chosen from {'$\texttt{Ball}$', '$\texttt{Box}$', '$\texttt{Halfspace}$', '$\texttt{Hyperplane}$'}. Each set involves corresponding parameters. For example, a box constraint requires specifying lower bound $l$ and upper bound $u$: when $l=-\infty$ and  $u=+\infty$, the box constraint reduces to an unconstrained case; when  $=0$ and $u=+\infty$, the box constraint becomes a nonnegativity constraint. For other types of $\Omega$, the relevant parameters are described in the model introduction above. Once $\Omega$ is chosen, $\texttt{SNSCO}$ can be invoked to solve the problem. The following Matlab code demonstrates the process of solving the recovery problem. 
</div>

<p style="line-height: 1;"></p>

```ruby
% demon recovery problems 
clc; close all; clear all; addpath(genpath(pwd));

K     = 10; 
M     = 10; 
N     = 100;
alpha = 0.05;
s     = ceil(alpha*N);

test = 2;  % Omega = {x|norm(x) <= r}  if test = 1
           % Omega = [lb,ub]^n         if test = 2    
           % Omega = {x|a'*x <= b}     if test = 3 
           % Omega = {x|Ax = b}        if test = 4 
sets = {'Ball','Box','Halfspace','Hyperplane'};
switch sets{test}
    case 'Ball'
        input1  = 2;
        input2  = [];
        xopt    = randn(K,1);
        xopt    = input1/norm(xopt)*xopt;
    case 'Box'
        input1  = -2; % can be -inf
        input2  = 2;  % can be inf
        xopt    = max(-10,input1)+(min(10,input2)-max(-10,input1))*rand(K,1); 
    case 'Halfspace'
        xopt    = rand(K,1);
        input1  = randn(K,1);
        input2  = sum(input1.*xopt)+rand(1); 
    case 'Hyperplane'
        xopt    = randn(K,1);
        input1  = randn(ceil(0.5*K),K);
        input2  = input1*xopt; 
end

% Generate data and define f and G
B      = randn(ceil(0.25*K),K)/sqrt(K);
d      = B*xopt;
BtB    = B'*B;
xi     = randn(K,M,N);
T      = randperm(N,s);
Mat    = rand(M,N);
D      = (Mat>=0.5) .* rand(M,N);
D(:,T) = (Mat(:,T)<1/3).*rand(M,nnz(T))-(Mat(:,T)>=2/3).*rand(M,nnz(T)); 
A      = reshape(xi,K,M*N);
C      = (squeeze(sum(xi .* xopt, 1))).^2 + D; 
Funcf  = @(x)FuncfRecovery(x,B,d,BtB);            % f(x)    = 0.5||Bx-d||^2
FuncG  = @(x,W,J)FuncGRecovery(x,W,J,A,C,K,M,N);  % G(x)_ij = <A_ij,x>^2-C_ij

% set parameters and call the solver
if  alpha  > 0.01
    pars.tau0 = 0.05+0.05*(test>4);
else
    pars.tau0 = 0.01;
    pars.thd  = 1e-1*(test==4)+1e-2*(test~=4);
end
out  = SNSCO(K,M,N,s,Funcf,FuncG,sets{test},input1,input2,pars);
fprintf(' Relative error: %7.3e \n', norm(out.x-xopt)/norm(xopt));
```

<div style="text-align:justify;">
The inputs and outputs of the Matlab version of $\texttt{SNSCO}$ are detailed below, with analogous specifications for the Python version. Inputs ($\texttt{K}$, $\texttt{M}$, $\texttt{N}$, $\texttt{s}$, $\texttt{Funcf}$, $\texttt{FuncG}$, $\texttt{FeaSet}$, $\texttt{input1}$, $\texttt{input2}$) are required. The parameters in  $\texttt{pars}$ are optional, but setting certain ones may improve the solver's performance and the solution quality.  It should be noted that $\texttt{FeaSet}$ can only be chosen from {'$\texttt{Ball}$', '$\texttt{Box}$', '$\texttt{Halfspace}$', '$\texttt{Hyperplane}$'}. For each set, the solver requires two inputs, $\texttt{input1}$ and $\texttt{input2}$. If an input is not needed, it can be set to empty $\texttt{[ ]}$. For example, when $\texttt{FeaSet}$ = '$\texttt{Ball}$', one may set $\texttt{input1}$ = 2 and $\texttt{input2}$ = $\texttt{[ ]}$ to represent a ball constraint with radius 2. When $\texttt{FeaSet}$ = '$\texttt{Box}$', one may set $\texttt{input1}$ = 0 and $\texttt{input2}$ = $\texttt{Inf}$ to represent a nonnegativity constraint.
</div>

<p style="line-height: 1;"></p>

```ruby
function out = SNSCO(K,M,N,s,Funcf,FuncG,FeaSet,input1,input2,pars)
% This solver solves 0/1 constrained optimization in the following form:
%
%         min_{x\in\R^K} f(x),  s.t. \| G(x) \|^+_0<=s, x\in Omega 
%
% where 
%      f(x) : \R^K --> \R
%      G(x) : \R^K --> \R^{M-by-N}
%      s << N 
%      \|Z\|^+_0 counts the number of columns with positive maximal values
%      Omega is a closed and convex set
% -------------------------------------------------------------------------
% Inputs:
%   K     : Dimnesion of variable x                              (REQUIRED)
%   M     : Row number of G(x)                                   (REQUIRED)
%   N     : Column number of G(x)                                (REQUIRED)
%   s     : An integer in [1,N), typical choice ceil(0.01*N)     (REQUIRED)
%   Funcf : Function handle of f(x)                              (REQUIRED)
%   FuncG : Function handle of G(x)                              (REQUIRED)
%   FeaSet: Feasible set for x, must be one of:                  (REQUIRED)
%          'Box'             [lb,ub]^K
%          'Ball'            {x|norm(x) <= r} 
%          'Halfspace'       {x|a'*x <= b}
%          'Hyperplane'      {x|Ax = b}
%           Default:         R^K
%   input1: A parameter related to FeasSet                       (REQUIRED)
%   input2: A parameter related to FeasSet                       (REQUIRED)
%   pars  : All parameters are OPTIONAL  
%           pars.x0    -- Initial point (default:ones(K,1)) 
%           pars.tau0  -- A vector of a number of  \tau0       (default  1)
%                         e.g.,pars.tau0=logspace(log10(.5),log10(1.75),20) 
%           pars.tol   -- Tolerance of halting condition (default 1e-6*M*N)
%           pars.maxit -- Maximum number of iterations       (default 2000) 
%           pars.disp  -- Show results or not at each step      (default 1)
% -------------------------------------------------------------------------
% Outputs:
%     out.x:      Solution x
%     out.obj:    Objective function value f(x)
%     out.G:      Function value of G(x) 
%     out.time:   CPU time
%     out.iter:   Number of iterations 
%     out.error:  Error
%     out.Error:  Error of every iteration
% -------------------------------------------------------------------------
% Written by Shenglong Zhou on 30/4/2024 based on the algorithm proposed in
%     Shenglong Zhou, Lili Pan, Naihua Xiu, and Geoffrey Ye Li, 
%     0/1 constrained optimization solving sample average approximation 
%     for chance constrained programming, Math Oper Res, 2024    	
% Send your comments and suggestions to <<< slzhou2021@163.com >>> 
% WARNING: Accuracy may not be guaranteed!!!!!  
% -------------------------------------------------------------------------
```
