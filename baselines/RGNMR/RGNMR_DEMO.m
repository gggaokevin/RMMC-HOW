clear all;
n1 = 500;
n2 = 500;
rank = 5;
condition_number = 5;
oversampling_ratio = 12;
alpha = 0.05;

% for non uniform sampling set this to true
power_sampling = false;

% generate L-star
singluar_values = linspace(1, 1/condition_number, rank) ;
[L_star, U, V, D, incoh] = generate_matrix(n1, n2, rank, singluar_values);

% generate omega
number_of_obsrved_entries = min(floor(rank*(n1+n2-rank)*oversampling_ratio), n1*n2);

if  ~power_sampling
    [mask, omega, flag]  = generate_mask(n1,n2 ,number_of_obsrved_entries, rank, 10000);
else
    beta = 2.5;
    [mask, omega, flag]  = power_law_mask(n1, n2, beta, number_of_obsrved_entries, rank, 10000);
end 

if flag == 0
    disp('FAILED TO GENERATE MASK');
    return
end

% gnerate corruption matrix S_star
S_tilde = unifrnd(-max(abs(L_star(:))), max(abs(L_star(:))), n1, n2);
[S_mask, lambda_star, flag] = generate_outliers(mask, omega, alpha, rank);
if flag == 0
    S_star = S_tilde .* S_mask;
else
    disp('FAILED TO GENERATE CORRUPTION MATRIX');
    return
end

% create the observable matrix X
X = mask .*(L_star + S_star);

% display RGNMR progression
opts.verbose = 1;            
% use spectral initialization
opts.init_option = 0;
% number of iterations
opts.max_outer_iter = 50;      
opts.max_inner_iter = 2000;

% stopping criteria
opts.stop_relRes = 1e-15;                                 
opts.stop_relDiff = 1e-15; 

[L_hat, convergence_flag,  all_relRes, iterations__with_unchanged_lambda] = ...
    RGNMR(X, omega, rank, length(lambda_star), opts);
fprintf("RELATIVE ERROR: %d\n", norm(L_hat - L_star, 'fro') / norm(L_star, 'fro'))