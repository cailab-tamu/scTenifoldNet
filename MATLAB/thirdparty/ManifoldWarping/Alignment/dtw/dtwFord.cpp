#include "mex.h"

#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define argmin(a, b)  (((a) < (b)) ? 0 : 1)

/* 
 * function [v, S, DC] = dtwFord(D)
 *
 * History
 *   create  -  Feng Zhou, 03-20-2009
 *   modify  -  Feng Zhou, 09-03-2010
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // D
    double *D = mxGetPr(prhs[0]);
    int n1 = mxGetM(prhs[0]);
    int n2 = mxGetN(prhs[0]);
    
    // v
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *v = mxGetPr(plhs[0]);

    // S
    plhs[1] = mxCreateDoubleMatrix(n1, n2, mxREAL);
    double *S = mxGetPr(plhs[1]);

    // DC
    plhs[2] = mxCreateDoubleMatrix(n1, n2, mxREAL);
    double *DC = mxGetPr(plhs[2]);

    // dynamic programming
    int dS[2][2] = {{1, 3}, {2, 3}};
    double d, tmpDC;
    int p, tmpS1, tmpS2;
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {    
            p = j * n1 + i;
            d = D[p];

            if (i == 0 && j == 0) {
                DC[p] = d;
                S[p] = 0;
            } else if (i == 0) {
                DC[p] = DC[p - n1] + d;
                S[p] = 1;
            } else if (j == 0) {
                DC[p] = DC[p - 1] + d;
                S[p] = 2;
            } else {
                tmpDC = min(DC[p - n1], DC[p - 1]);
                tmpS1 = argmin(DC[p - n1], DC[p - 1]);

                DC[p] = min(tmpDC, DC[p - n1 - 1]) + d;
                tmpS2 = argmin(tmpDC, DC[p - n1 - 1]);
                S[p] = dS[tmpS1][tmpS2];
            }
        }
    }
    *v = DC[n1 * n2 - 1];
}
