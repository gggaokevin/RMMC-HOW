function [M_out, U, V, RMSE_2, RMSE_2_out, iter, time] = L0_IQR(origin_X, X, Omega_array, rank, maxiter1, maxiter2, lambda, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RMC-l0 baseline = Sheng2025HOW Algorithm 1 ("HOW-IQR") with the Welsch-based
% outlier update replaced by the hard l0-proximal update of Liu2023 (the paper's
% own "l0-prox" variant: "replace (16) with (13), leaving the rest unchanged").
%
% Input
% -origin_X: ground truth (ONLY used when opts.gtStop = true, for diagnostics)
% -X: incomplete observation (zeros on missing entries)
% -Omega_array: sampling matrix (1 = observed)
% -rank: rank of origin_X
% -maxiter1/maxiter2: I1 / I2 = maximum iterations of the two stages
%                     (Sheng2025HOW convergence study: I1 = I2 = 100)
% -lambda: proximal parameter (Liu2023 released code: 1e-5)
% -opts (optional struct):
%     .xi1     threshold constant xi1 in  c_k = min(xi1*p_k, c_{k-1})
%              (Sheng2025HOW: xi1 = 2; p_k = IQR/1.349 = 0.7413*(q75-q25))
%     .c0      initial threshold c_0.  Sheng2025HOW Table I (p. 6, "Initial
%              threshold of HOW (c_0)") gives c_0 = 2, which is the value used
%              with eq. (28) c_{k+1} = min(xi1*p_{k+1}, c_k) -- hence the default
%              here is 2.  (In our normalized data xi1*p_1 << 2, so the cap does
%              not bind; opts.c0 = Inf reproduces the V5-T64 variant.)
%              The companion "initial kernel size" sigma_0 = 2*sqrt(2) of Table I
%              is unused here: the l0 variant replaces (16) by the hard-threshold
%              operator (13), so no kernel size enters.
%     .tol     termination tolerance on the relative reconstructed error
%              rel_e^k = (e^k - e^{k-1})/e^{k-1},  eps1 = eps2 = 1e-4
%              (Sheng2025HOW convergence study: "termination tolerance 1e-4")
%     .stopMode 'relchg' (default, ground-truth free) | 'gt' (legacy Liu code
%              stop: masked absolute RMSE to origin_X < 1e-4) | 'none' (run the
%              full I1/I2, as frozen in 实验方案_V5-T63 §0)
%     .stop1   whether stage 1 also tests the criterion (Algorithm 1 step 4,
%              default true; the legacy Liu code never stopped inside stage 1)
%
% e^k = || Theta_Omega - (X^k Y^k - Z^k) - O^k ||_F^2   (observed entries only;
% Theta_Omega is the INCOMPLETE data matrix, so the criterion never touches the
% ground truth -- Sheng2025HOW Sec. III text after eq. (30)).
%
% V5-T64 (2026): (a) c_0 = +Inf so that the threshold IS the signed IQR value
% xi1*0.7413*(q75-q25) from the first iteration on; (b) ground-truth-free
% termination per Algorithm 1 step 4 / step 6; (c) cap 100/100 by the caller.
% V5-A3 (2026): signed IQR convention (Robust2023 / RMMC) kept.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% -M_out   : M = U*V
% -U, -V   : factors
% -RMSE_2  : history of e^k (phase 2), zero-padded to maxiter2
% -RMSE_2_out : e^k at the last phase-2 iteration
% -iter    : TOTAL number of factor updates performed (phase 1 + phase 2)
% -time    : wall clock seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 8, opts = struct(); end
xi1      = getfield_def(opts,'xi1',2);
c        = getfield_def(opts,'c0',2);       % Sheng2025HOW Table I: c_0 = 2
tol      = getfield_def(opts,'tol',1e-4);    % V5-T64: eps1 = eps2 = 1e-4
stopMode = getfield_def(opts,'stopMode','relchg');
stop1    = getfield_def(opts,'stop1',true);

tic;
[m,n] = size(X); % row col
U = randn(m,rank);
V = randn(rank,n);
Z = randn(m,n);
S = Omega_array * 0;
RMSE_2 = zeros(1,maxiter2);
e_prev = NaN;
it1 = 0;

% ---------------- stage 1: plain factorization (O = 0) ----------------
for iter = 1 : maxiter1
    G = X + Z - S;
    U = (G * V' - lambda * U) * pinv(V * V' - lambda * eye(rank));
    V = pinv(U' * U - lambda * eye(rank)) * (U' * G - lambda * V);
    Z = U * V - U * V.*Omega_array;
    M = U * V;
    it1 = iter;

    e = norm((X - (M - Z) - S).*Omega_array,'fro')^2;
    if stop1 && stop_one(e_prev, e, tol, stopMode, origin_X, M-Z, Omega_array, m, n), break; end
    e_prev = e;
end

% ---------------- threshold c_1 from the signed IQR of the residual ----------
D = X + Z - M;
D_m_n = D(find(D));
d = 0.7413*(quantile(D_m_n,0.75)-quantile(D_m_n,0.25));   % = IQR/1.349
c = min([c xi1*d]);                                       % eq. (28): min(xi1*p_k, c_{k-1})
S = 0.*(abs(D)<=c)+ D.*(abs(D)>c);                        % hard l0 (Liu2023 (13))

% ---------------- stage 2: robust iterations with the l0-proximal O-update --
e_prev = NaN;
it2 = 0;
for iter = 1 : maxiter2
    G = X + Z - S;
    U = (G * V' - lambda * U) * pinv(V * V' - lambda * eye(rank));
    V = pinv(U' * U - lambda * eye(rank)) * (U' * G - lambda * V);
    Z = U * V - U * V.*Omega_array;
    M = U * V;
    it2 = iter;

    D = X + Z - M;
    D_m_n = D(find(D));
    d = 0.7413*(quantile(D_m_n,0.75)-quantile(D_m_n,0.25));
    c = min([c xi1*d]);

    S = 0.*(abs(D)<=c)+ D.*(abs(D)>c);

    e = norm((X - (M - Z) - S).*Omega_array,'fro')^2;
    RMSE_2(iter) = sqrt(e/(m*n));
    if stop_one(e_prev, e, tol, stopMode, origin_X, M-Z, Omega_array, m, n), break; end
    e_prev = e;
end
M_out = M;
RMSE_2_out = RMSE_2(it2);
iter = it1 + it2;      % total factor updates
time = toc;
end

% ---------------------------------------------------------------------------
function stop = stop_one(e_prev, e, tol, mode, origin_X, MZ, Om, m, n)
switch mode
    case 'gt'      % legacy Liu2023 released-code stop (Uses the ground truth)
        stop = norm((origin_X - MZ).*Om,'fro') / sqrt(m*n) < tol;
    case 'none'
        stop = false;
    otherwise      % 'relchg': ground-truth free, Sheng2025HOW Algorithm 1
        if isnan(e_prev) || e_prev <= 0
            stop = false;
        else
            stop = abs(e - e_prev)/abs(e_prev) < tol;
        end
end
end

function v = getfield_def(s, f, d)
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
