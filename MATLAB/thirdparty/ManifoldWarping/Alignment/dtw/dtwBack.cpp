#include "mex.h"
#include "stdio.h"

#define ptr(i, j, k) ((j) * k + i)

/* 
 * function [P, W] = dtwBack(S)
 *
 * History
 *   create  -  Feng Zhou, 03-20-2009
 *   modify  -  Feng Zhou, 08-04-2009
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // S
    double *S = mxGetPr(prhs[0]);
    int n1 = mxGetM(prhs[0]);
    int n2 = mxGetN(prhs[0]);

    // determine #step
    int p = ptr(n1 - 1, n2 - 1, n1);
    int n0 = 1;
    while (S[p] > 0) {
        if (S[p] == 1) {
            p -= n1;
        } else if (S[p] == 2) {
            --p;
        } else if (S[p] == 3) {
            p -= n1 + 1;
        } else {
            printf("unknown path direction\n");
            exit(1);
        }
        ++n0;
    }

    // P
    plhs[0] = mxCreateDoubleMatrix(n0, 2, mxREAL);
    double *P = mxGetPr(plhs[0]);

    // W
    plhs[1] = mxCreateDoubleMatrix(n1, n2, mxREAL);
    double *W = mxGetPr(plhs[1]);
    for (int i = 0; i < n1 * n2; ++i) {
        W[i] = 0;
    }

    // track back
    int i = n1 - 1, j = n2 - 1;
    p = ptr(i, j, n1);
    int n0a = n0;
    int n0b = n0 * 2;
    P[--n0a] = i + 1;
    P[--n0b] = j + 1;
    W[p] = 1;
    while (n0a > 0) {
        if (S[p] == 1) {
            --j;
            p -= n1;
            
        } else if (S[p] == 2) {
            --i;
            --p;
        } else if (S[p] == 3) {
            --i;
            --j;
            p -= n1 + 1;
            
        } else {
            printf("unknown path direction\n");
            exit(1);
        }
        
        P[--n0a] = i + 1;
        P[--n0b] = j + 1;
        W[p] = 1;
    }
}
