#include "mex.h"

//#define min(a, b)  (((a) < (b)) ? (a) : (b))
//#define argmin(a, b)  (((a) < (b)) ? 0 : 1)

/* 
 * function B = rowBd(P)
 *
 * History
 *   create  -  Feng Zhou, 03-20-2009
 *   modify  -  Feng Zhou, 09-03-2010
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // P
    double *P = mxGetPr(prhs[0]);
    int n0 = mxGetM(prhs[0]);
    int n1 = int(P[n0 - 1]);
    
    // B
    plhs[0] = mxCreateDoubleMatrix(n1, 2, mxREAL);
    double *B = mxGetPr(plhs[0]);

    int head = 0;
    int i, tail;
    while (head < n0) {
        i = int(P[head]) - 1;
        
        tail = head;
        while (tail < n0 && P[tail] - 1 == i) {
            ++tail;
        }
        
        B[i] = P[head + n0];
        B[i + n1] = P[tail - 1 + n0];
        head = tail;
    }
}
