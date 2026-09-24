function [convergence_flag] = check_early_convergence(iterations__with_unchanged_lambda, stop_lambda_convergence, all_relRes, stop_relRes, ...
                                                L_hat, L_hat_previous,  stop_relDiff, ...
                                                stop_relResDiff, iter, LSQR_iters_done, verbose)
    L_hat_diff = norm(L_hat - L_hat_previous, 'fro') / norm(L_hat, 'fro');
    convergence_flag = false; 
    if stop_lambda_convergence && iterations__with_unchanged_lambda > 1
            msg = '[INSIDE RGNMR] Early stopping: Lambda converged\n';
            convergence_flag = 1;
    elseif all_relRes(iter) < stop_relRes
        msg = '[INSIDE RGNMR] Early stopping: small error on observed entries\n';
        convergence_flag = 1;    
    elseif L_hat_diff < stop_relDiff
        msg = '[INSIDE RGNMR] Early stopping: L_hat does not change\n';
        convergence_flag = 1;
    elseif iter > 1 && ...
            abs(all_relRes(iter-1)/all_relRes(iter)-1) < stop_relResDiff
        msg = '[INSIDE RGNMR] Early stopping: relRes does not change\n';
        convergence_flag = 1;   
    elseif iter > 1 && LSQR_iters_done == 0
        msg = '[INSIDE RGNMR] Early stopping: no iterations of LSQR solver\n';
        convergence_flag = 1;
    end

    if convergence_flag
        if verbose
            fprintf(msg);
        end
    end
end