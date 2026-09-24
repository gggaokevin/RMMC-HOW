This directory contains the matlab implementation of RGNMR and demo files to demonstrate usage
# Directory Tree
+ utils
  + utils_RGNMR : functions used by RGNMR
    + utils_c : functions implemented in c
      + generate_sparse_matrix_indices.c
      + generate_sparse_matrix_indices.mexw64
      + generate_val_A.c
      + generate_val_A.mexw64
    + binary_weights.m
    + check_early_convergence.m
    + estimate_number_of_outliers.m
    + init_RGNMR.m
    + set_defaults_options.m
    + solve_LSQR_problem.m
    + vectorize_observed_matrix.m
  + test-utils : functions to generate synthetic data
    + generate_mask.m
    + generate_matrix.m
    + generate_outliers.m
    + power_law_mask.m
    
+ RGNMR.m
+ RGNMR_DEMO.m
+ RGNMR_BS_DEMO.m
