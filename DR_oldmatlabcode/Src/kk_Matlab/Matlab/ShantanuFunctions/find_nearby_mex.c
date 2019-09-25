/*=============================================================================
%FIND_NEARBY_MEX Find elements in target vector that are nearby each element of a reference vector.
%
%   IDX = FIND_NEARBY_MEX(X,Y,LOWER_BOUND,UPPER_BOUND) takes two column vectors
%   of real non-NaN double elements X and Y, and two double real scalars
%   LOWER_BOUND and UPPER_BOUND, and returns a column cell array IDX containing
%   index vectors. IDX is the same size as X, such that
%
%   IDX{i} = find((Y >= X(i) - LOWER_BOUND) & (Y < X(i) + UPPER_BOUND));
%   all(Y(IDX{i}) - X(i) < UPPER_BOUND)
%   all(Y(IDX{i}) - X(i) >= LOWER_BOUND)
%
%   This MEX function depends on standard C libraries.
%
%Written by SMK, 2009 May 29.
%
=============================================================================*/

#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  unsigned long i, j;

  double* x;
  size_t nx;
  double* y;
  size_t ny;
  double lower_bound;  
  double upper_bound;
  double nan_value;

  nan_value = mxGetNaN();

  /* Process arguments */
  if ( (nrhs != 4) || (nlhs > 1) ) {
    mexErrMsgTxt(
        "Usage: nearby_idx = find_nearby_mex(x,y,lower_bound,upper_bound)\n");
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
    mexErrMsgTxt("Third argument (lower bound) must be a non-NaN real "
        "double scalar\n");
  } else {
    lower_bound = mxGetScalar(prhs[2]);
  }
  if ( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
      (mxGetM(prhs[3]) != 1) || (mxGetN(prhs[3]) != 1) || 
      (mxGetNumberOfDimensions(prhs[3]) > 2) ||
      (mxGetScalar(prhs[3]) == nan_value) ) {
    mexErrMsgTxt("Fourth argument (upper bound) must be a non-NaN real "
        "double scalar\n");
  } else {
    upper_bound = mxGetScalar(prhs[3]);
  }
  if (!(upper_bound > lower_bound)) {
    mexErrMsgTxt("Fourth argument (upper bound) must be greater than "
        "third argument (lower bound)");
  }

  plhs[0] = mxCreateCellMatrix(nx, 1); 
  /* indices is a jagged array of mxArray objects, which will be used to
  populate plhs[0] */
  mxArray** indices;
  indices = (mxArray**) mxCalloc((mwSize) nx, (mwSize) sizeof(mxArray*));

  /* tmparray is an array of doubles */
  double* tmparray;
  /* Allocate tmparray to accomodate the maximum range of indices (at most,
  every element of Y could be listed in tmparray) */
  tmparray = (double*) mxCalloc((mwSize) ny, (mwSize) sizeof(double));

  /* count is the number of elements of tmparray that are to be copied to the
  output */
  unsigned long int count;

  for (i = 0; i < nx; i++) {
    /* Reset count */
    count = 0;
    /* Find elements of y that fall within [lower_bound, upper_bound) of x[i] */
    for (j = 0; j < ny; j++) {
      if ( (y[j] - x[i] >= lower_bound) && (y[j] - x[i] < upper_bound) ) {
        /* Add each matching element to the list of indices, using MATLAB
        convention of starting at 1 */
        tmparray[count++] = (double) (j + 1);
      }
    }
    /* Allocate indices[i] to accomodate count doubles */
    indices[i] = mxCreateDoubleMatrix((mwSize) count, 1, mxREAL);
    /* Copy the first count elements of tmparray to Pr* of indices[i] */
    memcpy((double*) mxGetPr(indices[i]), tmparray, 
        (size_t) (count*sizeof(double)));
    /* Assign indices[i] to be the ith cell of plhs[0] */
    mxSetCell(plhs[0], i, indices[i]);
  }

  /* Clean up */
  /* Don't destroy indices because plhs[0] points to it */
  mxFree(tmparray);
}

