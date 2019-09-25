/*=============================================================================
%FIND_NEARBY_MEX Find nearby values within a vector.
%
%   IDX = FIND_NEARBY_MEX(X,LOWER_BOUND,UPPER_BOUND) takes a column vector of
%   real non-NaN double elements X, a double real scalar LOWER_BOUND, and a
%   double real scalar UPPER_BOUND, and returns a column cell array IDX. IDX is
%   the same size as X, such that
%
%   IDX{i} = find((X >= X(i) - LOWER_BOUND) & (X < X(i) + UPPER_BOUND));
%
%   This MEX function depends on standard C libraries.
%
%Written by SMK, 2009 May 29.
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

  double* data;
  size_t numdata;
  double lower_bound;  
  double upper_bound;
  double nan_value;
  nan_value = mxGetNaN();

  /* Process arguments */
  if ( (nrhs != 3) || (nlhs > 1) ) {
    mexErrMsgTxt(
        "Usage: nearby_idx = find_nearby_mex(x,lower_bound,upper_bound)\n");
  }
  if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || 
      (mxGetNumberOfDimensions(prhs[0]) > 2) || (mxGetN(prhs[0]) != 1) ||
      (mxGetM(prhs[0]) < 1) ) {
    mexErrMsgTxt("First argument (data) must be a non-empty real double "
        "column vector\n");
  } else {
    data = mxGetPr(prhs[0]);
    numdata = mxGetNumberOfElements(prhs[0]);
  }
  for (i = 0; i < numdata; i++) {
    if (data[i] == nan_value) {
      mexErrMsgTxt("First argument (data) must not contain NaN values\n");
    }
  }  
  if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
      (mxGetM(prhs[1]) != 1) || (mxGetN(prhs[1]) != 1) || 
      (mxGetNumberOfDimensions(prhs[1]) > 2) || 
      (mxGetScalar(prhs[1]) == nan_value) ) {
    mexErrMsgTxt("Second argument (lower bound) must be a non-NaN real "
        "double scalar\n");
  } else {
    lower_bound = mxGetScalar(prhs[1]);
  }
  if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      (mxGetM(prhs[2]) != 1) || (mxGetN(prhs[2]) != 1) || 
      (mxGetNumberOfDimensions(prhs[2]) > 2) ||
      (mxGetScalar(prhs[2]) == nan_value) ) {
    mexErrMsgTxt("Third argument (upper bound) must be a non-NaN real "
        "double scalar\n");
  } else {
    upper_bound = mxGetScalar(prhs[2]);
  }
  if (!(upper_bound > lower_bound)) {
    mexErrMsgTxt("Third argument (upper bound) must be greater than "
        "second argument (lower bound)");
  }

  plhs[0] = mxCreateCellMatrix(numdata, 1); 
  /* indices is a jagged array of mxArray objects */
  mxArray** indices;
  indices = (mxArray**) mxMalloc((mwSize) (numdata*sizeof(mxArray*)));
  /* tmparray is an array of doubles */
  double* tmparray;
  /* Allocate tmparray to accomodate the maximum range of indices */
  tmparray = (double*) mxMalloc((mwSize) numdata*sizeof(double));
  /* count is the number of elements of tmparray that are to be copied to the
  output */
  unsigned long int count;

  for (i = 0; i < numdata; i++) {
    /* Reset count */
    count = 0;
    /* Find elements of data that fall within [lower_bound, upper_bound) of
    data[i] */
    for (j = 0; j < numdata; j++) {
      if ( (data[j] - data[i] >= lower_bound) && 
          (data[j] - data[i] < upper_bound) ) {
        /* Add to list of indices, using MATLAB convention of starting at 1 */
        tmparray[count++] = j + 1;
      }
    }
    indices[i] = mxCreateDoubleMatrix((mwSize) count, 1, mxREAL);
    /* Copy the first count elements of tmparray to the real data part of
    indices[i] */
    memcpy((double*) mxGetPr(indices[i]), tmparray, 
        (size_t) (count*sizeof(double)));
    /* Assign indices[i] to be the ith cell of plhs[0] */
    mxSetCell(plhs[0], i, indices[i]);
  }

  /* Clean up */
  mxFree(tmparray);
}

