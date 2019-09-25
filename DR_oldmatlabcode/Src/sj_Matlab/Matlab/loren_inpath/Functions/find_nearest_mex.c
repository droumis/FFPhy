/*=============================================================================
%FIND_NEAREST_MEX Find element in target vector closest to specified offset relative to each element of reference vector.
%
%   [IDX, DELTA] = FIND_NEAREST_MEX(X,Y,OFFSET) takes two column vectors of real
%   non-NaN double elements X and Y, and a double real scalar OFFSET, and
%   returns column vectors IDX and DELTA of the same size as X, such that
%
%   DELTA = X + OFFSET - Y(IDX)
%   abs(DELTA(i)) = min(abs(X(i) + OFFSET - Y)) for all i = 1:numel(X)
%
%   If the minimum is not unique, then IDX(i) is chosen to be as close to i as
%   possible (which is useful behavior when X is identical to Y).
%
%   This MEX function depends on standard C libraries.
%
%Written by SMK, 2009 November 25.
%
=============================================================================*/

#define _FILE_OFFSET_BITS 64
#define fopen fopen64
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  mwSize alloc_sz;
  unsigned long i, j;

  double* x;
  size_t nx;
  double* y;
  size_t ny;
  double offset;
  double nan_value;
  nan_value = mxGetNaN();
  double inf_value;
  inf_value = mxGetInf();

  /* Process arguments */
  if ( (nrhs != 3) || (nlhs > 2) ) {
    mexErrMsgTxt(
        "Usage: [idx, delta] = find_nearest_mex(x,y,offset)\n");
  }
  if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || 
      (mxGetNumberOfDimensions(prhs[0]) > 2) || (mxGetN(prhs[0]) != 1) ||
      (mxGetM(prhs[0]) < 1) ) {
    mexErrMsgTxt("First argument (X) must be a non-empty real double "
        "column vector\n");
  } else {
    x = mxGetPr(prhs[0]);
    nx = mxGetNumberOfElements(prhs[0]);
  }
  for (i = 0; i < nx; i++) {
    if (x[i] == nan_value) {
      mexErrMsgTxt("First argument (X) contains NaN\n");
    }
  }  
  if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || 
      (mxGetNumberOfDimensions(prhs[1]) > 2) || (mxGetN(prhs[1]) != 1) ||
      (mxGetM(prhs[1]) < 1) ) {
    mexErrMsgTxt("Second argument (Y) must be a non-empty real double "
        "column vector\n");
  } else {
    y = mxGetPr(prhs[1]);
    ny = mxGetNumberOfElements(prhs[1]);
  }
  for (j = 0; j < ny; j++) {
    if (y[j] == nan_value) {
      mexErrMsgTxt("Second argument (Y) contains NaN\n");
    }
  }  
  if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      (mxGetM(prhs[2]) != 1) || (mxGetN(prhs[2]) != 1) || 
      (mxGetNumberOfDimensions(prhs[2]) > 2) || 
      (mxGetScalar(prhs[2]) == nan_value) ) {
    mexErrMsgTxt("Third argument (target offset) must be a non-NaN real "
        "double scalar\n");
  } else {
    offset = mxGetScalar(prhs[2]);
  }

  double* indices;
  double* deltas;
  plhs[0] = mxCreateDoubleMatrix(nx, 1, mxREAL); 
  indices = mxGetPr(plhs[0]);
  if (nlhs == 2) {
    plhs[1] = mxCreateDoubleMatrix(nx, 1, mxREAL);
    deltas = mxGetPr(plhs[1]);
  } else {
    deltas = mxMalloc((mwSize) nx*sizeof(double));
  }
  /* Initialize indices and delta with Inf values so that comparisons
  automatically fail */
  for (i = 0; i < nx; i++) {
    indices[i] = inf_value;
    deltas[i] = inf_value;
  }

  double tmp_delta;

  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) {
      tmp_delta = y[j] - (x[i] + offset);
      if ( fabs(tmp_delta) < fabs(deltas[i]) ) {
        deltas[i] = tmp_delta;
        /* Convert C-style index starting from zero to MATLAB-style index
        starting from 1 */
        indices[i] = (double) (j + 1);
      } else if ( fabs(tmp_delta) == fabs(deltas[i]) ) {
        if ( fabs(((double) j) - ((double) i)) < 
            fabs(indices[i] - 1 - ((double) i)) ) {
          deltas[i] = tmp_delta;
          /* Convert C-style index starting from zero to MATLAB-style index
          starting from 1 */
          indices[i] = (double) (j + 1);
        }
      }
    }
  }

  /* Clean up if we mxMalloc'ed deltas apart from plhs[1] */
  if (nlhs != 2) {
    mxFree(deltas);
  }

}

