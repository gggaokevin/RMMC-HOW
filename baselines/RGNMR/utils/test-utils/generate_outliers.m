function [H, lambda_star_2d, reject] = generate_outliers(mask, omega, alpha, rank)
    %% Input
    % mask - n1xn2 sparse matrix 
    %     ones at chosen  entries and
    %     zeros at unchosen  entries.
    % omega - all indices of chosen enries in mask 
    % alpha - fraction of outliers 

    %% Output
    reject = 1;
    reject_number = 0;
    reject_max = 10;
    n1 = size(mask,1);
    n2 = size(mask, 2);
    while(reject == 1)
        lambda_star = (sort(randperm(size(omega, 1), floor(alpha*size(omega, 1)))))';
        lambda_star_2d = omega(lambda_star,:);
        i_lambda_star = lambda_star_2d(:,1);
        j_lambda_star = lambda_star_2d(:,2);
        H = sparse(i_lambda_star,j_lambda_star,ones(floor(alpha*size(omega, 1)),1),n1,n2);
        nr_entr_col_W_H = sum(mask-H,1)';
        nr_entr_row_W_H = sum(mask-H,2);

        if (isempty(find(nr_entr_col_W_H<rank,1)) == 1) && (isempty(find(nr_entr_row_W_H<rank,1)) == 1)
            reject = 0;
        else
            reject_number = reject_number + 1;
            if reject_number > reject_max
                 disp('No mask found!');
                 break
            end 
        end 
    end
end