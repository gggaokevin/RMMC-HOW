function [vectorized_X] = vectorize_observed_matrix(X, omega)    
    number_of_observed_entries = size(omega,1);   % number of observed entries
    vectorized_X = zeros(number_of_observed_entries,1);
    for counter=1:number_of_observed_entries
            vectorized_X(counter) = X(omega(counter,1),omega(counter,2));
    end
end
