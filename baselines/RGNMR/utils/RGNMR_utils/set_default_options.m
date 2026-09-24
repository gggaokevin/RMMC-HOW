function [opts] = set_default_options(opts)
%% configurations: default values for option variables
opts_default.verbose = 0;                       % display intermediate results

% number of iterations and tolerance
opts_default.max_outer_iter = 100;              % maximal number of outer iterations
opts_default.max_inner_iter = 2000;             % maximal number of inner iterations for the LSQR solver
opts_default.inner_init_tol = 1e-15;             % initial tolerance of the LSQR solver
opts_default.LSQR_smart_tol = 0;                % use LSQR_tol==relRes*1e-1 

% initialization
opts_default.init_option = 0;                   % 0 for SVD, 1 for random, 2 for opts.init_U, opts.init_V
opts_default.init_U = NaN;                      % if opts.init_option==2, use this initialization for U
opts_default.init_V = NaN;                      % if opts.init_option==2, use this initialization for V

% early stopping criteria (-1 to disable a criterion)
opts_default.stop_relRes = 1e-16;               % small relRes threshold (relevant to noise-free case)
opts_default.stop_relDiff = -1;                 % small relative X_hat difference threshold
opts_default.stop_relResDiff = -1;              % small relRes difference threshold
opts_default.stop_relResStuck_ratio = -1;       % stop if minimal relRes didn't change by a factor of stop_relResNoChange_ratio...
opts_default.stop_relResStuck_iters = -1;       % ... in the last #stop_relResNoChange_iters outer iterations

% additional configurations for the binary search of the number of outliers
opts_default.estimates_outliers = 0;            % when estimating the amount of outliers set to true (1), otherwise set to false (0) 
opts_default.failing_criterion = 10e-4;         % the required threshold to consider the amount of outliers removed as sufficent 
opts_default.stop_lambda_convergence = false;   % if true the algorithm stops if the estimate of the set of outliers converged


% for each unset option set its default value
opts_default.stop_D_diff = -1;   % V5 patch: field used by RGNMR.m but missing in release
fn = fieldnames(opts_default);
for k=1:numel(fn)
    if ~isfield(opts,fn{k}) || isempty(opts.(fn{k}))
        opts.(fn{k}) = opts_default.(fn{k});
        fn{k} = opts_default.(fn{k});
    end
end
end

