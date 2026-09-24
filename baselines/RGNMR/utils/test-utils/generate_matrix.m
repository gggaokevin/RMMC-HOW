function [L_star, U, V, D, mu] = generate_matrix(n1, n2, rank, singluar_values)

%% Input 
% n1 - number of rows
% n2 - number of columns
% rank - the rank of the matrix toi generate
% singular_values - the singular values of the matrix

%% Output
% L_star = UDV it's SVD decomposition
% mu : the incohernece of the matrix 


D = diag(singluar_values);


Z = randn(n1,rank); 
[U, ~, ~] = svd(Z,'econ'); 

Z = randn(n2,rank); 
[V, ~, ~] = svd(Z,'econ'); 

L_star = U * D * V';

% find the incohernece 
mu_U = max(sum(U.^2,2)) * n1 / rank;
mu_V = max(sum(V.^2,2)) * n2 / rank;

mu = max(mu_U, mu_V);

end