#include "mex.h"
#include <string.h>

void generate_sparse_A(double *U, double *V, double *omega, mwSize nv, mwSize r, mwSize n1, mwSize n2, double *colind_A, double *val_A) {
    mwSize i, j, k;
    mwSize idx;
    
    for (i = 0; i < nv; i++) {
        j = (mwSize)omega[i] - 1;
        k = (mwSize)omega[i + nv] - 1;
        
        for (idx = 0; idx < r; idx++) {
            val_A[i * 2 * r + idx] = U[j + idx * n1];
            val_A[i * 2 * r + r + idx] = V[k + idx * n2];
        }
    }
}

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *U, *V, *omega, *colind_A, *val_A;
    mwSize nv, r, n1, n2;
    
    /* Check for proper number of arguments */
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:generate_sparse_A:invalidNumInputs", "Four inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MATLAB:generate_sparse_A:invalidNumOutputs", "One output required.");
    }
    
    /* Get the inputs */
    U = mxGetPr(prhs[0]);
    V = mxGetPr(prhs[1]);
    omega = mxGetPr(prhs[2]);
    colind_A = mxGetPr(prhs[3]);
    
    /* Get the dimensions of the inputs */
    nv = mxGetM(prhs[2]);
    r = mxGetN(prhs[0]);
    n1 = mxGetM(prhs[0]);
    n2 = mxGetM(prhs[1]);
    
    /* Create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1, 2 * r * nv, mxREAL);
    val_A = mxGetPr(plhs[0]);
    
    /* Call the computational routine */
    generate_sparse_A(U, V, omega, nv, r, n1, n2, colind_A, val_A);
}
