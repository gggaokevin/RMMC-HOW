function [L_hat, convergence_flag,  all_relRes, iterations_since_D_changed] = ...
    RGNMR(X,omega, rank, outliers_num, opts, varargin)

%
% Written by  Eilon Vaknin Laufer & Boaz Nadler / 2025
% based on the code of GNMR by Pini Zilber and Boaz Nadler  
% with modifications by Christian Kuemmerle (sparse matrix construction of A)
%
% INPUT: 
% X = the observed matrix
% omega = list of pairs (i,j) of the observed entries
% rank = target rank of the underlying matrix
% outliers_num = the expected number of outlires
% opts = options meta-variable (see set_default_options for details)
%
% OUTPUT:
% L_hat = approximation of the underlying matrix
% all_relRes = list of the residual errors of the least-squares problem throughout the iterations
% iter = final iteration number
% convergence_flag = indicating whether the algorithm converged
%

%% set RGNMR configuration
opts  = set_default_options(opts);

%% initalize RGNMR
% for a detailed explnation about the matrix D see solve_LSQR_problem
[U, V, L_hat, D] = init_RGNMR(opts.init_option, opts.init_U, opts.init_V, X, omega, rank, outliers_num);

%% define variables for later use 
[~, n2] = size(X);   % number of rows and colums
colind_A = generate_sparse_matrix_indices(omega, rank, n2); % used when constructing the sparse matrix A

%% set variables before iterations
all_relRes = zeros(opts.max_outer_iter,1); % stores the relRes of all iterations
convergence_flag = 0; % if true breaks from iteration 
number_of_iterations_since_D_changed = 0;
LSQR_tol = opts.inner_init_tol; % the tolorence used by the LSQR solver

%% iterations
for iter = 1:opts.max_outer_iter
    
    L_hat_previous = L_hat;
    D_previous = D;
  
    %% solve the LSQR problem
    [U, V, L_hat, entriwise_residulas, relRes, LSQR_iters_done] = solve_LSQR_problem(X, U, V,omega, D, colind_A, LSQR_tol, opts.max_inner_iter);
    all_relRes(iter) = relRes; % Record before invoking the release stopping helper.

    %%  update the estimate of the set of non corrupted entries represented by D
    D = binary_weights(entriwise_residulas,  outliers_num);
    number_of_iterations_since_D_changed  = (number_of_iterations_since_D_changed + 1) *(isequal(D, D_previous));

    
    %% decrease the tolorence
    if opts.LSQR_smart_tol
        LSQR_tol = max(2*eps, relRes*1e-1);
    end
    
    %% report
    if opts.verbose
        X_hat_diff = norm(L_hat - L_hat_previous, 'fro') / norm(L_hat, 'fro');
        fprintf('[INSIDE GNMR] iter %4d \t diff X_r %5d\t relRes %6d\n',...
            iter, X_hat_diff, relRes);
    end
    %% check early convergence
    % Use the published MATLAB helper's argument order and configured criteria.
    % D is the binary weighting matrix for the estimated uncorrupted set;
    % the helper calls its unchanged-set counter "unchanged_lambda".
    convergence_flag = check_early_convergence( ...
        number_of_iterations_since_D_changed, opts.stop_lambda_convergence, ...
        all_relRes, opts.stop_relRes, L_hat, L_hat_previous, ...
        opts.stop_relDiff, opts.stop_relResDiff, iter, LSQR_iters_done, opts.verbose);
    if convergence_flag
        break
    end   
end

iterations_since_D_changed = number_of_iterations_since_D_changed;
%% return
[U_hat, lambda_hat, V_hat] = svds(L_hat ,rank);
L_hat = U_hat * lambda_hat * V_hat';
end
