
#include <mex.h>
#include <stdint.h>
#include <algorithm>

class arg_cmp {  // comparator for argsort
public:
    double* row;
    arg_cmp(double* r): row(r) {}
    bool operator()(int a, int b){ return row[a]<row[b]; }
};

void knnsearch_same(const mxArray* X, const int k, int32_t* idx) {
    const int m = mxGetM(X);
    const int n = mxGetN(X);
    // D = repmat(sum(X.*X)',[1 n]) + repmat(sum(X.*X),[n 1]) - 2*X'*X;
    double** D = new double*[n];
    for (int i=0; i<n; ++i) {
        D[i] = new double[n];
        D[i][i] = 0;
    }
    if (mxIsSparse(X)) {
        mwIndex* jc = mxGetJc(X);
        mwIndex* ir = mxGetIr(X);
        double* pr = mxGetPr(X);
        for (int i=0; i<n; ++i) {
            for (int j=i+1; j<n; ++j) {
                D[i][j] = 0;
                // only look at nonzero elements in columns i and j
                mwIndex p1 = jc[i];
                mwIndex p2 = jc[j];
                const mwIndex p1max = jc[i+1];
                const mwIndex p2max = jc[j+1];
                while (p1 < p1max || p2 < p2max) {
                    double x1 = 0;
                    double x2 = 0;
                    if (p2 == p2max || (p1 < p1max && ir[p1] <= ir[p2]))
                        x1 = pr[p1++];
                    if (p1 == p1max || (p2 < p2max && ir[p1] >= ir[p2]))
                        x2 = pr[p2++];
                    const double diff = x1 - x2;
                    D[i][j] += diff*diff;
                }
                D[j][i] = D[i][j];
            }
        }
    } else {
        // full matrices are so nice and straightforward
        double* x = mxGetPr(X);
        for (int i=0; i<n; ++i) {
            for (int j=i+1; j<n; ++j) {
                D[i][j] = 0;
                for (int l=0; l<m; ++l) {
                    const double diff = x[i*m+l] - x[j*m+l];
                    D[i][j] += diff*diff;
                }
                D[j][i] = D[i][j];
            }
        }
    }
    // [~,idx] = sort(D); idx = idx(2:k+1,:)';
    int* tmp_idx = new int[n];
    for (int i=0; i<n; ++i) {
        // always reorder because partial_sort mutates the array
        for (int j=0; j<n; ++j)
            tmp_idx[j] = j;
        // skip the first one, it's always 0 for D[i][i]
        std::partial_sort(tmp_idx, tmp_idx+k+1, tmp_idx+n, arg_cmp(D[i]));
        for (int j=1; j<k+1; ++j) {
            idx[(j-1)*n+i] = tmp_idx[j] + 1;  // silly matlab, with your 1-based indexes
        }
        delete[] D[i];
    }
    delete[] tmp_idx;
    delete[] D;
}

void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){
    // function idx = knnsearch_same(X,k)
    if (nlhs != 1)
        mexErrMsgTxt("Must assign result to something");
    if (nrhs != 2)
        mexErrMsgTxt("Invalid args.\nUsage: idx = knnsearch_same(X,k)");
    
    const int n = mxGetN(prhs[0]);
    const int k = (int) mxGetScalar(prhs[1]);
    if (k >= n)
        mexErrMsgTxt("Only n-1 nearest neighbors exist");
    plhs[0] = mxCreateNumericMatrix(n, k, mxINT32_CLASS, mxREAL);
    int32_t* idx = (int32_t*)mxGetData(plhs[0]);
    knnsearch_same(prhs[0], k, idx);
}
