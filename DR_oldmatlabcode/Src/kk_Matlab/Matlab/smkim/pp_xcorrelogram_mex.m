%PP_XCORRELOGRAM_MEX Compute the raw cross-correlogram (joint event histogram) of a bivariate point process realization.
%
%   CC = PP_XCORRRELOGRAM_MEX(TARGET,REFERENCE,LAGS,BIN_WIDTH) takes the event
%   times TARGET and REFERENCE and computes the raw (uncorrected)
%   cross-correlogram of TARGET events in bins of size BIN_WIDTH centered at
%   LAGS. All arguments must be MATLAB double class. TARGET, REFERENCE and LAGS
%   must be vectors, and BIN_WIDTH must be a positive scalar.
%
%   This MEX function depends on standard C libraries. It is meant to be called
%   through the wrapper M-function PP_XCORRELOGRAM.
%
%
%Written by SMK, 2009 October 1.
%
