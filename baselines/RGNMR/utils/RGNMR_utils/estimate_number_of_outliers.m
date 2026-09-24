function [num_of_outliers] = estimate_number_of_outliers...
(X, omega, rank, alpha_min, alpha_max, opts, test)

%% INPUT
% X - observed matrix 
% omega - list of pairs (i,j) of the observed entries
% rank - target rank of the underlying matrix
% alpha_min - lower bound on the fraction of corrupted entries
% alpha_max - upper bound on the fraction of corrupted entries
% opts - options meta-variable (see set_deafult_options method for details)
% test - optinal, for testing 

%% OUTPUT
% num_of_outliers - upper bound on the number of outliers

if nargin < 8
    test = false;
end

%% set options for RGNMR
opts.stop_lambda_convergence = true;
j = 1;
alpha = alpha_min;

%% Binary Search fro the fraction of outliers 
while((alpha_max - alpha_min)*length(omega) > 1)
    [~,  ~, ~, iterations_with_unchanged_outliers_set] = RGNMR(X, omega, rank, ceil(alpha*length(omega)), opts, test);
     
     % if Lambda converged 
     if iterations_with_unchanged_outliers_set > 0
         alpha_min = alpha;
         alpha = (alpha_max + alpha) / 2;
     else % if Lambda did not converge 
         alpha_max = alpha;
         alpha = (alpha_min + alpha) / 2;
     end
     j = j+1;
     if opts.verbose_num_of_outliers
        fprintf("[ESTIMATING NUMBER OF OUTLIERS]k_min  = %d,    k_max = %d\n", ceil(alpha_min*length(omega)), ceil(alpha_max*length(omega)))
     end
end

alpha = alpha_max;
num_of_outliers = ceil(alpha*length(omega));
end