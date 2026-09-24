function [Out_X, S, U, V, info] = HOW_RMC(M, Omega_array, rak, maxiter, xi1, xi2)
% HOW_RMC -- Robust Matrix Completion via Hybrid Ordinary-Welsch (RMC-HOW).
% Faithful port of the OFFICIAL released code of
%   Z.-Y. Wang, H. C. So, A. M. Zoubir, "Robust low-rank matrix recovery via
%   hybrid ordinary-Welsch function," IEEE TSP 2023  (RMC_HOW.m, SASD version)
% from RMMC_reference/RMC-HOW-HOC-HOP-main/.
% -------------------------------------------------------------------------
% This is a SINGLE-MATRIX (per-frame) robust MC solver. For a multi-frame
% sequence apply it frame-by-frame (see HOW_RMC_multi below / callers).
%
% Inputs:
%   M          : m x n observed matrix (missing entries = 0)
%   Omega_array: m x n logical/0-1 sampling mask (1 = observed)
%   rak        : rank of the object matrix
%   maxiter    : max outer iterations (robust phase)
%   xi1, xi2   : HOW threshold multipliers: c = xi1*loc, sigma = xi2*loc.
%                (Robust release default: xi1=2, xi2=sqrt(2)*xi1.)
% Outputs:
%   Out_X : recovered low-rank matrix (U*V)
%   S     : per-entry outlier component (sparse)
%   U,V   : factors
%   info  : struct with loc/c/sigma trajectory (diagnostics)
% -------------------------------------------------------------------------
if nargin < 5, xi1 = 2;  end
if nargin < 6, xi2 = xi1*sqrt(2); end
if nargin < 4, maxiter = 100; end
if islogical(Omega_array), Omega_array = double(Omega_array); end

Pro  = @(x,c,scale) 0.*(abs(x)<=c) + (abs(x)-abs(x).*exp((c^2-x.^2)./scale^2.*(abs(x)>c))).*sign(x).*(abs(x)>c);
f_welsch = @(x,c,sigma) 0.5*x.^2.*(abs(x)<=c) + (sigma^2/2*(1-exp((c^2-x.^2)/sigma^2))+c^2/2).*(abs(x)>c);

[m,n] = size(M);
S = 0;
loc = 1000;
U = randn(m,rak);
V = randn(rak,n);
X = U*V;

% ---------- phase 1: plain factorization warm start (official: up to 100) ----
warmup = 0;
for k = 1:100
    warmup = k;
    D = M - S;
    dU = -(D - Omega_array.*(U*V))*V';
    du = -dU*((V*V')^(-1));
    tu = -trace(dU'*du)/(norm(Omega_array.*(du*V),'fro'))^2;
    U = U + tu*du;
    dV = -U'*((D - Omega_array.*(U*V)));
    dv = (U'*U)^(-1)*dV;
    tv = -trace(dV'*dv)/(norm(Omega_array.*(U*dv),'fro'))^2;
    V = V + tv*dv;
    X_new = U*V;
    X_err = norm(X_new-X,'fro')^2/norm(X,'fro')^2;
    X = X_new;
    if X_err < 1e-3, break; end
end

% ---------- phase 2: robust HOW iterations (outlier filtering + SASD) --------
loc_traj = zeros(1,maxiter); c_traj = zeros(1,maxiter); s_traj = zeros(1,maxiter);
loss_HOW = [];
for iter = 1:maxiter
    T = M - X.*Omega_array;                 % observed residual
    tv_ = T(T~=0);
    loc_1 = 0.7413*(quantile(tv_,0.75) - quantile(tv_,0.25));   % = IQR/1.349
    loc = min(loc, loc_1);                  % monotone shrink from 1000
    scale = xi2*loc;
    c = xi1*loc;
    S = Pro(T, c, scale);                   % outlier support estimate
    loc_traj(iter)=loc; c_traj(iter)=c; s_traj(iter)=scale;

    D = M - S;
    dU = -(D - Omega_array.*(U*V))*V';
    du = -dU*((V*V')^(-1));
    tu = -trace(dU'*du)/(norm(Omega_array.*(du*V),'fro'))^2;
    U = U + tu*du;
    dV = -U'*((D - Omega_array.*(U*V)));
    dv = (U'*U)^(-1)*dV;
    tv = -trace(dV'*dv)/(norm(Omega_array.*(U*dv),'fro'))^2;
    V = V + tv*dv;
    X = U*V;
    loss_HOW(end+1) = sum(sum(f_welsch(M - X.*Omega_array, c, scale))); %#ok<AGROW>
    if numel(loss_HOW) >= 2 && abs(loss_HOW(end)-loss_HOW(end-1))/abs(loss_HOW(end-1)) < 1e-4
        break;                              % official early stop (loss rel. change)
    end
end
Out_X = X;
info = struct('loc',loc_traj(1:iter),'c',c_traj(1:iter),'sigma',s_traj(1:iter),'iter',iter, ...
              'warmup',warmup);   % V5-T70: stage-1 count so that total iters = warmup + iter
end
