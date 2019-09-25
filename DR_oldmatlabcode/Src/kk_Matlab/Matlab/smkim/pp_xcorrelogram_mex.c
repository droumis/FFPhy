/*=============================================================================
%PP_XCORRELOGRAM_MEX Compute the raw cross-correlogram (joint event histogram) of a bivariate point process realization.
%
%   CC = PP_XCORRRELOGRAM_MEX(REFERENCE,TARGET,LAGS,BIN_WIDTH) takes the event
%   times REFERENCE and TARGET and computes the raw (uncorrected)
%   cross-correlogram of TARGET events in bins of size BIN_WIDTH centered at
%   LAGS. All arguments must be MATLAB double class. REFERENCE, TARGET and LAGS
%   must be vectors, and BIN_WIDTH must be a positive scalar.
%
%   This MEX function depends on standard C libraries. It is meant to be called
%   through the wrapper M-function PP_XCORRELOGRAM.
%
%
%Written by SMK, 2009 October 1.
%
=============================================================================*/


#include "mex.h"

void usage();

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  /* Check input arguments */
  if ( !(nrhs == 4) || (nlhs > 1) ) {
    usage();
  }
  if ( !mxIsNumeric(prhs[0]) || (mxGetNumberOfDimensions(prhs[0]) != 2) || 
      !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || 
      !((mxGetM(prhs[0]) == 1) || (mxGetN(prhs[0]) == 1)) ) {
    mexErrMsgTxt("REFERENCE must be a vector of real doubles\n");
  }
  if ( !mxIsNumeric(prhs[1]) || (mxGetNumberOfDimensions(prhs[1]) != 2) || 
      !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || 
      !((mxGetM(prhs[1]) == 1) || (mxGetN(prhs[1]) == 1)) ) {
    mexErrMsgTxt("TARGET must be a vector of real doubles\n");
  }
  if ( !mxIsNumeric(prhs[2]) || (mxGetNumberOfDimensions(prhs[2]) != 2) || 
      !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) || mxIsEmpty(prhs[2]) ||
      !((mxGetM(prhs[2]) == 1) || (mxGetN(prhs[2]) == 1)) ) {
    mexErrMsgTxt("LAGS must be a non-empty vector of real doubles\n");
  }
  if ( !mxIsNumeric(prhs[3]) || (mxGetNumberOfDimensions(prhs[3]) != 2) || 
      !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) || mxIsEmpty(prhs[3]) ||
      (mxGetNumberOfElements(prhs[3]) != 1) || !(mxGetScalar(prhs[3]) > 0) || 
      mxIsNaN(mxGetScalar(prhs[3])) || mxIsInf(mxGetScalar(prhs[3])) ) {
    mexErrMsgTxt("BIN_WIDTH must be a positive finite scalar\n");
  }

  double* reference = mxGetPr(prhs[0]);
  size_t n_reference = mxGetNumberOfElements(prhs[0]);
  double* target = mxGetPr(prhs[1]);
  size_t n_target = mxGetNumberOfElements(prhs[1]);
  double* lags = mxGetPr(prhs[2]);
  size_t n_lags = mxGetNumberOfElements(prhs[2]);
  double bin_width = mxGetScalar(prhs[3]);
  /* mxCreateDoubleMatrix initializes the data block to all zeros */
  plhs[0] = mxCreateDoubleMatrix(1, (mwSize) n_lags, mxREAL);
  double* cc = mxGetPr(plhs[0]);

  unsigned long i, j, k;
  double lag_ij;
  for (i = 0; i < n_reference; i++) {
    for (j = 0; j < n_target; j++) {
      lag_ij = target[j] - reference[i];
      for (k = 0; k < n_lags; k++) {
        if ((lag_ij >= lags[k] - bin_width/2) && 
            (lag_ij < lags[k] + bin_width/2)) {
          cc[k]++;
        }
      }
    }
  }

} 

void usage() {
  mexErrMsgTxt(
"   CC = PP_XCORRELOGRAM_MEX(REFERENCE,TARGET,LAGS,BIN_WIDTH)\n"
"     REFERENCE, TARGET, LAGS must be double vectors\n"
"     BIN_WIDTH must be positive double scalar\n"
  );
}


