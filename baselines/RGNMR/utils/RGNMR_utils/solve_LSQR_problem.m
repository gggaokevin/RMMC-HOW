function [U, V, L_hat, entriwise_residuals, relRes, LSQR_iters_done] = ...
         solve_LSQR_problem(X, U, V,omega, Lambda, colind_A, tol, max_iter)

    %% INPUT:
    % X - observed matrix
    % U, V - current factor matrices 
    % omega - list of pairs (i,j) of the observed entries
    % Lambda - projection matrix to the estimated set of non corrupted entries
    % colind_A - column indices of non zero entries in A
    % tol - tolorence for the lsqr MATLAB method
    % max_iter - maximal number of iteration for the lsqr MATLB method

    %% OUTPUT:
    % U, V - updated factor matrices
    % L_hat - current estimate of L*
    % entriwise_residuals - (L_hat - X) as a vector projected to 
    %                       the set of observed entries
    % relRes - relres output of the lsqr MATLAB method
    % LSQR_iters_done -  iteration number iter at which lsqr terminated

    %% Explanation
    % We construct the least of squares(lsqr) problem in eq.1 (see paper). 
    % as a weighted lsqr problem. 
    % The diagonal weights matrix Lambda has 0 or 1 on the diagonal. 
    % 0 corresponds to entries suspected to be outliers. 
    % We then solve min ||Lambda*(Ax-b)||_F.
    % We construct the updated U, V from the solution x.

    
    [n1, ~] = size(U);
    [n2, r] = size(V);
    %% construct variables A,b for LSQR to solve |DAx-Db|^2 
    % construct A
    A = generate_sparse_A(U, V, omega, colind_A);    
    % construct b 
    X_updated = X + U*V';
    b = vectorize_observed_matrix(X_updated, omega);

    %% solve the least squares problem
    [x, ~, relRes, LSQR_iters_done] = lsqr(Lambda*A, Lambda*b, tol, max_iter);

    %% construct U_next and V_next from the solution x
    U_next = zeros(size(U));
    V_next = zeros(size(V)); 
    nc_list = r * (0:1:(n2-1)); 
    for i=1:r
        V_next(:,i) = x(i+nc_list); 
    end
    nr_list = r * (0:1:(n1-1)); 
    start_idx = r*n2; 
    for i=1:r
        U_next(:,i) = x(start_idx + i + nr_list);
    end
    L_hat = U * V_next' + U_next * V' - U * V';
    U = U_next;
    V = V_next;
    entriwise_residuals = abs(A*x - b);
end

%% auxiliary function for generating sparse matrix A

function A = generate_sparse_A(U,V,omega,colind_A)
    
    %% INPUT:
    % U, V - current factor matrices
    % omega -  list of pairs (i,j) of the observed entries
    % colind_A - column indices of non zero values in A
    
    %% OUTPUT:
    % A  - sparse matrix used to construct the lsqr problem

    
    nv = size(omega,1); % number of observed entries
    rank = size(U,2);
    n1 = size(U,1);
    n2 = size(V,1);
    
    % generate row indices of non zero entires in A
    rowind_A=kron(1:nv,ones(1,2*rank));

    % extract non zero values in A
    val_A = generate_val_A(U, V, omega, colind_A);

    % map val_A into (rowind_A, colind_A) in an nv X (rank*(n1+n2)) matrix
    A = sparse(rowind_A, colind_A, val_A, nv, rank*(n1+n2));
end
