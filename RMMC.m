%% Robust Multi-matrix Completion (RMMC) -- V6 / snapshot "RMMC_V561"
%% Frozen 2026-09-23 against main.tex of the same date: Sec. 3.1-3.2 and
%% Algorithm 1 (eqs. (1)-(9)).  Everything reported so far was produced by an
%% earlier fingerprint; this file is the frozen reference for the next round.
function [USV, O, Err_1, i, Snap] = RMMC(IncompleteData, array_Omega_c, CompleteData, r, c, max_out_iter, eta1, eta2, lambda, mode, stop_tol, snap_at, rho1_in, rho2_in, verbose)
% mode (optional, default 0): 0 = HOW proximal (the paper's model);
% 1 = ell_0 hard threshold at the same threshold; 2 = ell_1 soft threshold.
% Modes 1 and 2 are DIAGNOSTIC ABLATIONS ONLY and are not used for any reported
% result (every experiment in the paper runs mode 0).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RMMC: MMC (Sheng et al., TPAMI 2026) + HOW robust channel (Wang-So-Zoubir,
%%% Robust low-rank matrix recovery via hybrid ordinary-Welsch, IEEE TSP 2023).
%%%
%%% THIS FILE IS THE V6 IMPLEMENTATION, aligned with Algorithm 1 and Sec. 3.1-3.2
%%% of main.tex.  Notation mapping paper -> code:
%%%     M_k  -> Y        N_k  -> Z        G_k -> G        W_k -> Q
%%%     O_k  -> O        U,V  -> U,V      S_k -> S        A_k -> A
%%%     X_Omega^c -> X_Omega_c            eta_1,eta_2 -> eta1,eta2
%%% (The multipliers were called P,Q in the V5 code; main.tex now denotes them
%%%  G,W.  The threshold multipliers were called xi1,xi2; main.tex uses eta_1,eta_2.
%%%  Argument positions are unchanged, so every existing caller keeps working.)
%%%
%%% Objective (paper eq. (3)) and augmented Lagrangian (paper eq. (4)):
%%%   min sum_k [ 1/2 ||X~_k - U S_k V^T - O_k||_{Omega_k}^2 + Phi(O_k)
%%%               + lambda ||L U S_k V^T||_F^2 + lambda ||U S_k V^T R||_F^2 ]
%%% Every block update below is the exact minimiser of (4); the data term carries
%%% NO factor 2 (unlike the MMC baseline, whose data term is ||.||^2 and therefore
%%% keeps the 2(...) numerator of the released MMC.m).
%%%
%%% Round i of Algorithm 1, in this file's order:
%%%   [step 2] E_k^i = Pi_{Omega_k}(X~_k - M_k^{i-1})                 [eq. (5)]
%%%            d_k^i = IQR(vec(E_k^i))/1.349, evaluated on the OBSERVED
%%%            support only: the mask-induced zeros of Omega_k^c are NOT in the
%%%            quantile sample (observed-support convention, as in the released
%%%            HOW code of Wang2023HOW);
%%%            O_k^i = p_{c_k^{i-1},s_k^{i-1}}(E_k^i)  -- PREVIOUS thresholds;
%%%            c_k^i = min(c_k^{i-1}, eta1 d_k^i), s_k^i = min(s_k^{i-1}, eta2 d_k^i).
%%%            c_k^0 = s_k^0 = +inf, hence O_k^1 = 0: the robust channel opens
%%%            in round 2, as Sec. 3.1 states.  The +inf sentinels are never used in
%%%            an arithmetic expression: round 1 is handled by an explicit flag, so a
%%%            degenerate scale d_k^1 = 0 cannot produce eta*inf = NaN.
%%%   [step 3] A_k^1 = X~_k and A_k^i = M_k^{i-1} + G_k^{i-1}/rho1 for i >= 2.
%%%   [step 4] three warm-started GLRAM sweeps on A_k^i give U^i, V^i and
%%%            S_k^i = U^i' A_k^i V^i.
%%%   [step 5] the closed-form blocks of eq. (6), followed by the multiplier
%%%            updates of eq. (6):
%%%              M_k^i = H_L^{-1}[ X_k^{i-1} - O_k^i - G_k^{i-1} + W_k^{i-1}
%%%                                  + rho1 U^i S_k^i V^i' + rho2 N_k^{i-1} ],
%%%              X_k^i   = X~_k + Pi_{Omega_k^c} M_k^i,
%%%              N_k^i = (rho2 M_k^i - W_k^{i-1}) H_R^{-1},
%%%              G_k^i = G_k^{i-1} + rho1 (M_k^i - U^i S_k^i V^i'),
%%%              W_k^i = W_k^{i-1} + rho2 (N_k^i - M_k^i),
%%%            with H_L = (1+rho1+rho2) I + 2 lambda L'L and
%%%                 H_R = rho2 I + 2 lambda R R'.
%%%            There is NO observed-position recalibration of G and no re-pinning
%%%            of A: the iteration is exactly Algorithm 1.
%%%
%%% Convergence-side identities these updates satisfy (Sec. 3.2), useful when
%%% auditing the code against the proof:
%%%     W_k^i = -2 lambda N_k^i R R'                                  (exact)
%%%     G_k^i = Pi_{Omega_k}(X~_k - O_k^i) + Pi_{Omega_k^c} M_k^{i-1}
%%%               - K M_k^i + N_k^{i-1}(rho2 I - 2 lambda R R'),
%%%     with K = (1+rho2) I + 2 lambda L'L.
%%%
%%% Degenerate boundary (safety net, not expected): if at least three quarters of the
%%% OBSERVED residuals of Omega_k were exactly zero, both quartiles would vanish
%%% and d_k^i = 0.  Sec. 3.1 assumes d_k^i > 0; the first round would then fix
%%% c = s = 0 and the HOW expression of the next round would evaluate 0/0.  The run is
%%% therefore stopped at that point and the cell is flagged through a NEGATIVE returned
%%% iteration index (i < 0), with USV/O/Err_1 left at the rounds actually completed.
%%% No silent repair of the scale is applied.
%%%
%%% These update rules fix the iterates, so every reported number must come from
%%% this file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  IncompleteData: 2D/3D incomplete data
%         array_Omega_c: I - sampling_matrix (1 = missing)
%         CompleteData: groundtruth data for relative error
%         r: row dimension of reduced representation S
%         c: column dimension of reduced representation S
%         max_out_iter: max number of iterations, 100 is suggested
%         eta1, eta2: HOW threshold multipliers (eta_1, eta_2 of the paper); FIXED
%                   in 03_code (the selection of the (eta1, eta2) pair is the
%                   subject of the selection study).  Historical callers pass them
%                   under the names xi1, xi2; the argument positions are identical.
%         lambda (optional): Tikhonov smoothing weight; the final value is selected
%                   in 03_code/s4_hp_select.m
% Output: USV: reconstruction
%         O:   per-frame outlier matrices of the LAST completed round
%         Err_1: relative reconstruction error per outer iteration
%         i: current outer iteration (frame count is K, frame index is k);
%            i < 0 means the run was stopped because the IQR interval collapsed to
%            zero (mask-inclusive scale, missing ratio >= 75%), i.e. the cell is not
%            executable at this setting.
%         Snap: snapshots requested through snap_at
%-------------------------------------------------------------------
[m,n,K] = size(IncompleteData);            % needed by the data-family default rule
[lam_def, rho1_def, rho2_def] = rmmc_defaults(m,n);
% Defaults below are the REPORTED settings for the data family of this input,
% switched on the slice geometry by rmmc_defaults.m:
%   synthetic S1-S3 (100 x 200): lambda = 5e-3
%   ORL       S4    (112 x  92): lambda = 0.1
% both families share the reported pair (rho1,rho2) = (9.05,1.45)
% Explicit arguments always override the rule.
if nargin < 9 || isempty(lambda), lambda = lam_def; end
if nargin < 10 || isempty(mode), mode = 0; end
% stop_tol:
%   > 0 : GT-FREE stopping on the relative change of the ADMM iterate M,
%         relchg = mean_k ||M_k^i - M_k^{i-1}||_F / ||M_k^i||_F < stop_tol.
%         The DEFAULT is the paper's 1e-4.
%   = 0 : the LEGACY GT-BASED monitor (Err_1(i) < 1e-3) stays ACTIVE.  This is NOT
%         "early stop disabled": it is a GROUND-TRUTH stop (independent audit
%         V5-T46).  It must never decide a reported run; it is reachable only if a
%         caller asks for stop_tol = 0 explicitly.
if nargin < 11 || isempty(stop_tol), stop_tol = 1e-4; end
% snap_at (V5-T71, optional): list of outer-iteration indices whose reconstruction
% is to be returned in the 5th output Snap (struct array with fields iter, USV).
% Used to report RMMC@It100 alongside RMMC@It200 from ONE 200-iteration run:
% the iteration is deterministic (fixed init, no RNG in the loop), so the iterate
% captured at i=100 is bit-identical to running max_out_iter=100. Default [] keeps
% every existing caller bit-identical and adds no memory.
if nargin < 12 || isempty(snap_at), snap_at = []; end
% rho1_in / rho2_in (optional): the ADMM penalties of the two splitting
% constraints.  Defaults come from rmmc_defaults.m and follow the same data-family
% switch as lambda; both pairs are the smallest satisfying the sufficient
% conditions (eq. (9) of the paper) for the respective lambda.
if nargin < 13 || isempty(rho1_in), rho1_in = rho1_def; end
if nargin < 14 || isempty(rho2_in), rho2_in = rho2_def; end
% verbose (optional, default true): the per-iteration trace.  Default true keeps
% every pre-existing call identical; the hyper-parameter sweep turns it off because
% the parallel workers would otherwise flood the client console.
if nargin < 15 || isempty(verbose), verbose = true; end
rho1 = rho1_in;
rho2 = rho2_in;
% ([m,n,K] were resolved at the top, before the data-family default rule.)

% L: (m-1) x m and R: n x (n-1) first-difference operators (the "L" of eq. (3);
% the released MMC.m carries the same two lines under a stale "laplacian" comment)
L = [-eye(m-1), zeros(m-1,1)] + [zeros(m-1,1), eye(m-1)]; % L
R = [-eye(n-1); zeros(1,n-1)] + [zeros(1,n-1); eye(n-1)]; % R

% HOW proximal operator (paper eq. (2), real-valued data)
p_phi = @(D, c_, s_) sign(D) .* max(abs(D) - abs(D) .* exp((c_.^2 - D.^2) ./ (s_.^2)), 0);

% Initialisation (Algorithm 1): M_k = N_k = G_k = W_k = O_k = 0 and X_k = X~_k.
G   = zeros(m,n,K);   % G_k : multiplier of M_k = U S_k V^T
Q   = zeros(m,n,K);   % W_k : multiplier of N_k = M_k
Z   = zeros(m,n,K);   % N_k : column-smoothness copy
Y   = zeros(m,n,K);   % M_k : fidelity / row-smoothness copy (M_k^0 = 0)
M_prev = zeros(m,n,K);% M_k^{i-1}: input of the robust block of round i
X_Omega_c = zeros(m,n,K); % Pi_{Omega^c} X
A   = IncompleteData; % A_k^1 = X~_k (the i = 1 branch of step 3)
U   = [eye(r); zeros(m-r, r)]; % GLRAM initial L0
O   = zeros(m,n,K);   % outlier matrices (O_k^1 = 0)

% HOW thresholds: c_k^0 = s_k^0 = +inf, which gives O_k^1 = 0.  Round 1 is
% detected by its index, not by testing the sentinel, so no eta*inf product is ever
% formed.
c_thr = inf(1,K);
s_thr = inf(1,K);

% Preallocation
S = zeros(r,c,K);
USV = zeros(m,n,K);

GLRAM_ITER = 3; % GLRAM max iteration num

Err_1(1:max_out_iter) = zeros;
Snap = struct('iter', {}, 'USV', {}, 'elapsed', {}, 'setup', {});   % V5-T71/82: snapshots
t_start = tic;                                         % V5-T72/82: elapsed at snapshot, and setup

% Loop-invariant block inverses (depend only on m,n,rho,lambda,L,R):
% H_L = (1+rho1+rho2) I + 2 lambda L'L   and   H_R = rho2 I + 2 lambda R R'
inv_C_C = inv( (1 + rho1 + rho2) * eye(m) + 2 * lambda * L' * L );
inv_D_D = inv( rho2 * eye(n) + 2 * lambda * R * R');
setup_s = toc(t_start);   % V5-T82: RMMC's OWN one-off preparation (init + block inverses), reported explicitly

% Loop-invariant observed-support mask and per-frame ground-truth norm
% (relative-error denominator), computed once instead of every iteration.
Om = (array_Omega_c == 0);
gt_norm = zeros(1,K);
for k = 1:K, gt_norm(k) = norm(CompleteData(:,:,k),'fro'); end

for i = 1:max_out_iter % outer iteration

    if verbose, fprintf('RMMC Iteration: %d/%d\n', i, max_out_iter); end

    % ---------------- Robust block: step 2 of Algorithm 1 (eq. (5)) -----------
    % Run at the START of round i on M_k^{i-1} (M_k^0 = 0, so E_k^1 = Pi(X~_k)),
    % so that O_k^i is available to the M-update of the SAME round i.
    for k = 1:K
        res = (IncompleteData(:,:,k) - M_prev(:,:,k)) .* Om(:,:,k); % E_k^i
        % Robust scale on the vectorised masked residual (paper Sec. 3.1):
        %   d = IQR(vec(E))/1.349, with 0.7413 = 1/1.349 the unbiased sigma scale
        % for real-valued residuals.  The sample is res(Om(:,:,k)), i.e. the
        % OBSERVED positions only: the structural zeros that the mask puts on
        % Omega^c are excluded, which is the observed-support convention of the
        % released HOW code (Wang2023HOW) that this paper follows.
        res_k = res(Om(:,:,k));
        q1 = quantile(res_k, 0.25); q3 = quantile(res_k, 0.75);
        % Degenerate boundary (safety net, not expected): should at least three
        % quarters of the OBSERVED residuals be exactly zero, both quartiles vanish
        % and d_k^i = 0.  The first round then fixes c = s = 0,
        % and the HOW expression p_phi of the NEXT round would evaluate the
        % indeterminate form 0/0 in exp((c^2 - e^2)/s^2).  Paper Sec. 3.1 assumes
        % d_k^i > 0; instead of silently repairing the estimate we stop the run and
        % flag the cell, leaving USV/O/Err_1 at their values for the rounds completed
        % (i.e. at the initialisation when the failure happens in round 1).
        if q1 == 0 && q3 == 0
            i = -abs(i);   % negative iteration index = this cell is NOT executable
            if verbose
                fprintf(['RMMC: empty IQR interval at round %d (d = 0 over the ' ...
                         'observed support); cell skipped\n'], abs(i));
            end
            return;
        end
        d = 0.7413 * (q3 - q1);
        if i == 1
            % c_k^0 = s_k^0 = inf  =>  O_k^1 = 0 (robust channel opens in round 2)
            O(:,:,k) = zeros(m,n);
            c_thr(k) = eta1 * d;
            s_thr(k) = eta2 * d;
        else
            c_old = c_thr(k); s_old = s_thr(k);   % c_k^{i-1}, s_k^{i-1}
            O(:,:,k) = p_phi(res, c_old, s_old);  % p_phi(0)=0 on Omega^c automatically
            if mode == 1           % ell_0 ablation variant (same threshold convention)
                O(:,:,k) = res .* (abs(res) >= c_old);
            elseif mode == 2       % ell_1 ablation variant
                O(:,:,k) = sign(res) .* max(abs(res) - c_old, 0);
            end
            % Thresholds advance AFTER O, as in eq. (5): O_k^i uses c_k^{i-1}.
            c_thr(k) = min(eta1 * d, c_old);
            s_thr(k) = min(eta2 * d, s_old);
        end
    end
    % ------------------------------------------------------------------------

    % X~_k - O_k^i: the outlier-corrected observation of the M-block data term of
    % eq. (6).  (O is zero on Omega^c by construction, so Xc is zero there.)
    Xc = IncompleteData - O;

    % ---------------- Step 3: the GLRAM input A_k^i ------------------------
    if i >= 2
        A = G / rho1 + M_prev;   % A_k^i = M_k^{i-1} + G_k^{i-1}/rho1
    end
    % (i = 1 keeps A = X~_k, set at initialisation.)

    % ---------------- Step 4: factor block -----------------------------------
    % Three warm-started GLRAM sweeps on A_k^i, started from the previous U.
    [U,V]  = GLRAM(A, U, r, c, GLRAM_ITER, K);

    err(1:K) = zeros;
    if stop_tol > 0, M_ref = Y; end   % GT-free stop reference: M_k^{i-1}

    % ---------------- Step 5: remaining blocks and multipliers (eq. (6)) ------
    for k = 1:K % inner iteration over frames

        S(:,:,k) = U' * A(:,:,k) * V;

        USV(:,:,k) = U * S(:,:,k) * V';

        Y(:,:,k) = inv_C_C * ( ( Xc(:,:,k) + X_Omega_c(:,:,k))...
            - G(:,:,k) + rho1 * USV(:,:,k) + Q(:,:,k) + rho2 * Z(:,:,k) );

        G(:,:,k) = G(:,:,k) + rho1 * ( Y(:,:,k) - USV(:,:,k) );

        X_Omega_c(:,:,k) = Y(:,:,k) .* array_Omega_c(:,:,k);   % Pi_{Omega^c} X_k^i

        Z(:,:,k) = ( rho2 * Y(:,:,k) - Q(:,:,k)) * inv_D_D;

        Q(:,:,k) = Q(:,:,k) + rho2 * (Z(:,:,k) - Y(:,:,k));

        % relative reconstruction error
        err(k) = norm(CompleteData(:,:,k) - USV(:,:,k),'fro') / gt_norm(k);

    end

    if any(snap_at == i)   % V5-T71/72: snapshot + ACTUAL elapsed time at iteration i
        Snap(end+1) = struct('iter', i, 'USV', USV, 'elapsed', toc(t_start), 'setup', setup_s); %#ok<AGROW>
    end

    % M_k^i produced by this sweep is the input of the robust block of round i+1
    % (E_k^{i+1} = Pi_{Omega_k}(X~_k - M_k^i)), per eq. (5).
    M_prev = Y;

    Err_1(i) = sum(err)/K;

    if stop_tol > 0
        relchg = 0;
        for k = 1:K
            relchg = relchg + norm(Y(:,:,k)-M_ref(:,:,k),'fro')/max(norm(Y(:,:,k),'fro'),eps);
        end
        relchg = relchg/K;
        if relchg < stop_tol
            fprintf('RMMC stopped (GT-free) at iteration %d (relchg=%.2e, relerr=%.6f)\n', i, relchg, Err_1(i));
            break;
        end
    elseif Err_1(i) < 0.001 % 1e-3: LEGACY GT-BASED stop (ACTIVE when stop_tol==0; uses ground truth!)
        fprintf('RMMC converged at iteration %d\n', i)
        break;
    end

end
end
