function [Lambda, missed_outliers] = binary_weights(entriwise_residulas, num_of_outliers, omega, test)
    % V5 patch (2026): omega/test optional -- release main loop (RGNMR.m L50)
    % calls with 2 args while the signature declares 4; reference read-only.
    if nargin < 3, omega = []; end
    if nargin < 4, test = []; end
    %
    % binary_weights returns a diagonal matrix Lambda with ones and zeros on the diagonal.
    % Lambda(i,i) = 0 iff entriwise_residulas(i) is one of the num_of_outliers 
    % largest values in entriwise_residulas.
    %
    %% INPUT:
    %   entriwise_residulas - residuals vectors
    %   num_of_outliers: Number of largest rows to consider
    %
    %% OUTPUT:
    %   Lambda: Diagonal matrix with ones and zeros on the diagonal

    
    missed_outliers = 0;
    % Get the num_of_outliers indices with largest residuals 
    [C, indices] = sort(entriwise_residulas, 'descend');
    
    indices_of_outliers = indices(1:num_of_outliers);
    
    % if in test print number of missed outliers
    if isstruct(test)
        fprintf("missed outliers %d\n", length(setdiff(test.s, omega(indices_of_outliers, :),'rows')))
        missed_outliers = length(setdiff(test.s, omega(indices_of_outliers, :), 'rows'));
    end

    %Construct the diagonal matrix Lambda with ones and zeros
    Lambda = speye(length(entriwise_residulas));
    Lambda(indices_of_outliers, indices_of_outliers) = 0;
end