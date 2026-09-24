#include "mex.h"

void generate_sparse_matrix_indices(double *omega, int r, int n2, int nv, double *colind_A) {
    int counter, j, k, idx1, idx2;
    for (counter = 0; counter < nv; counter++) {
        j = (int)omega[counter];      // Cast to int
        k = (int)omega[nv + counter]; // Cast to int
        for (int i = 0; i < r; i++) {
            idx1 = r * (k - 1) + i;
            idx2 = r * n2 + r * (j - 1) + i;
            colind_A[counter * 2 * r + i] = idx1 + 1; // MATLAB is 1-based index
            colind_A[counter * 2 * r + r + i] = idx2 + 1; // MATLAB is 1-based index
        }
    }
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *omega; /* Input matrix */
    double *colind_A; /* Output matrix */
    int r; /* Input scalar */
    int n2; /* Input scalar */
    int nv; /* Number of elements in omega */
    
    /* Check for proper number of arguments */
    if(nrhs != 3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }
    
    /* Get the inputs */
    omega = mxGetPr(prhs[0]);
    r = (int) mxGetScalar(prhs[1]);
    n2 = (int) mxGetScalar(prhs[2]);
    nv = mxGetM(prhs[0]);
    
    /* Create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1, 2 * r * nv, mxREAL);
    colind_A = mxGetPr(plhs[0]);
    
    /* Call the function */
    generate_sparse_matrix_indices(omega, r, n2, nv, colind_A);
}
