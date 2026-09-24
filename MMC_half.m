%% MMC baseline in the SAME 1/2-normalized form as RMMC (internal baseline)
function [USV, Err_1, i] = MMC_half(IncompleteData, array_Omega_c, CompleteData, r, c, max_out_iter, lambda, stop_tol)
% V4-internal MMC BASELINE for fair comparison with RMMC.m. Identical to the
% released MMC (MMC_lambda.m, parent-paper ||.||^2 form) EXCEPT the data
% fidelity is written with the SAME 1/2 normalization as RMMC's HOW objective
% (eq:rmmc-o): the Y-step denominator is (1+rho1+rho2)I and the data term
% carries no factor 2. Consequently lambda here has the RMMC meaning
% (lambda_RMMC); run it at the SAME lambda as RMMC so that data/smoothing
% balance AND the ADMM penalty scaling (rho=1) coincide exactly for the two
% methods. (Parent-paper ||.||^2 protocol lambda = 2 * lambda_RMMC.)
% Kept in dev/: the released MMC.m / dev/MMC_lambda.m stay untouched as the
% parent-code reference.
% lambda (required, RMMC scale): 0.1 images / 0.15 video / 5e-4 synthetic.
if nargin < 7 || isempty(lambda), lambda = 0.1; end
% V5-T37: GT-free early stopping (same rule as RMMC.m): stop when the relative
% change of the ADMM iterate Y between outer iterations is below stop_tol.
if nargin < 8 || isempty(stop_tol), stop_tol = 0; end
rho1 = 0.5; rho2 = 0.5;   % V5-T51: rho = rho_MMC/2 (penalty-to-objective ratio kept as in MMC eq.(8))
[m,n,K] = size(IncompleteData);
L = [-eye(m-1), zeros(m-1,1)] + [zeros(m-1,1), eye(m-1)];
R = [-eye(n-1); zeros(1,n-1)] + [zeros(1,n-1); eye(n-1)];
P = zeros(m,n,K); Q = zeros(m,n,K); Z = zeros(m,n,K);
X_Omega_c = zeros(m,n,K);
A = IncompleteData;
U = [eye(r); zeros(m-r, r)];
S = zeros(r,c,K); USV = zeros(m,n,K); Y = zeros(m,n,K);
GLRAM_ITER = 3;
Err_1(1:max_out_iter) = zeros;
inv_C_C = inv( (1 + rho1 + rho2) * eye(m) + 2 * lambda * L' * L );
inv_D_D = inv( rho2 * eye(n) + 2 * lambda * R * R');
gt_norm = zeros(1,K);
for k = 1:K, gt_norm(k) = norm(CompleteData(:,:,k),'fro'); end
for i = 1:max_out_iter
    [U,V] = GLRAM(A, U, r, c, GLRAM_ITER, K);
    err(1:K) = zeros;
    if stop_tol > 0, Y_prev = Y; end
    for k = 1:K
        S(:,:,k) = U' * A(:,:,k) * V;
        USV(:,:,k) = U * S(:,:,k) * V';
        Y(:,:,k) = inv_C_C * ( ( IncompleteData(:,:,k) + X_Omega_c(:,:,k))...
            - P(:,:,k) + rho1 * USV(:,:,k) + Q(:,:,k) + rho2 * Z(:,:,k) );
        P(:,:,k) = P(:,:,k) + rho1 * ( Y(:,:,k) - USV(:,:,k) );
        X_Omega_c(:,:,k) = Y(:,:,k) .* array_Omega_c(:,:,k);
        Z(:,:,k) = ( rho2 * Y(:,:,k) - Q(:,:,k)) * inv_D_D;
        Q(:,:,k) = Q(:,:,k) + rho2 * (Z(:,:,k) - Y(:,:,k));
        A(:,:,k) = P(:,:,k) / rho1 + Y(:,:,k);
        A(:,:,k) = IncompleteData(:,:,k) + A(:,:,k) .* array_Omega_c(:,:,k);
        err(k) = norm(CompleteData(:,:,k) - USV(:,:,k),'fro') / gt_norm(k);
    end
    Err_1(i) = sum(err)/K;
    if stop_tol > 0
        relchg = 0;
        for k = 1:K
            relchg = relchg + norm(Y(:,:,k)-Y_prev(:,:,k),'fro')/max(norm(Y(:,:,k),'fro'),eps);
        end
        relchg = relchg/K; Y_prev = Y;
        if relchg < stop_tol
            fprintf('MMC_half stopped (GT-free) at iteration %d (relchg=%.2e, relerr=%.6f)\n', i, relchg, Err_1(i));
            break;
        end
    elseif Err_1(i) < 0.001
        fprintf('MMC_half converged at iteration %d\n', i)
        break;
    end
end
end
