function [omega, omega_2d, flag] = power_law_mask(n1, n2, beta, w, rank, max_resapmles)
    
    %% Generate power-law distributed values
    u = linspace(1, n1, n1);
    raw = (u).^(-1 / (beta - 1));  
    x = w * raw / sum(raw);  
    
    u = linspace(1, n2, n2);
    raw = (u).^(-1 / (beta - 1));  
    y = w * raw / sum(raw); 

    P = (x' * y) / w;


    reject_counter = 0;
    reject = true;  
    while reject && reject_counter < max_resapmles

        omega = sparse(rand(n1, n2) < P);
        nv = sum(omega(:));
        nr_entr_col = sum(omega,1)';
        nr_entr_row = sum(omega,2);
    
        if (isempty(find(nr_entr_row<rank+1,1)) == 0) || (isempty(find(nr_entr_col<rank+1,1)) == 0)
            reject_counter=reject_counter+1;
        else
            reject = false;
        end

        if reject_counter >= max_resapmles
            disp('No mask found!');
            break
        end
    end
    
    [idx_i, idx_j] = find(omega);
    omega_2d = [idx_i, idx_j];
    flag = ~reject;
end