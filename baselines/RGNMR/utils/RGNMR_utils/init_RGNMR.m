function [U, V, L_hat, Lambda] = init_RGNMR(init_option, U_init, V_init, X, omega, rank, outliers_num, test)
    % V5 patch (2026): official release calls with 7 args but body uses `test`;
    % guard nargin so the local copy runs (reference kept read-only).
    if nargin < 8, test = []; end
    %% INPUT:
    % init_option - if 0 apllay threshold operator then use spectral initilaization.
    %               if 1 use random initialization 
    %               if 2 use user-defined matrices U_init and V_init
    % X - observed matrix
    % omega - list of pairs (i,j) of the observed entries
    % rank - rank of the target matrix
    % outlier_num - an upper bound on the number of outliers in X
    % test - optinal, for testing

    %% OUTPUT:
    % U, V - initialization for factor matrices
    % L_hat - initial estimate, UV', of the target matrix L*
    % Lambda - projection matrix to the estimated set of non corrupted entries
    
    [n1, n2] = size(X);
    %% initialize U and V (of sizes n1 x r and n2 x r)
    if init_option == 0
        % initialization by rank-r SVD of observed matrix, after applying a
        % threshold operator (remove_top_fraction)
        [U, ~, V] = svds(remove_top_fraction(X,  2*outliers_num/(n1*n2)), rank);
    elseif init_option == 1
        % initialization by random orthogonal matrices
        Z = randn(n1,rank);
        [U, ~, ~] = svd(Z,'econ'); 
        Z = randn(n2,rank);
        [V, ~, ~] = svd(Z,'econ'); 
    else
        % initiazliation by user-defined matrices
        U = U_init;
        V = V_init; 
    end
    
    % compute intial estimate of L*
    L_hat = U * V';

    % vectorize the input matrix X
    vector_X = vectorize_observed_matrix(X, omega);
    % vectorize the inital estimator
    vector_L_hat = vectorize_observed_matrix(L_hat, omega);
    % construct an estimate of the set of non corrupted entries
    Lambda = binary_weights(abs(vector_L_hat - vector_X), outliers_num, omega, test);
end

function A = remove_top_fraction(A, alpha)

    %% Input:
    % A - a matrix of size n1 X n2
    % alpha - upper bound on the fraction of corrupted entries

    %% Output:
    % zeros the alpha*n2 largest entries in each row.
    % zeros the alpha*n1 largest entries in each column.

    [n1, n2] = size(A);
    num_row = max(1, round(alpha * n2));
    num_col = max(1, round(alpha * n1));

    % remove the largest magnitude entries in each row
    for i = 1:n1
        [~, idx] = sort(abs(A(i, :)), 'descend');
        A(i, idx(1:num_row)) = 0;
    end

    % remove the largest magnitude entries in each column
    for j = 1:n2
        [~, idx] = sort(abs(A(:, j)), 'descend');
        A(idx(1:num_col), j) = 0;
    end
end