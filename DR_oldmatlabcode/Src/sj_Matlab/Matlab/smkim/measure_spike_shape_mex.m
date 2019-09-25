%MEASURE_SPIKE_SHAPE Measure shape of spike waveforms recorded in NSpike .tt format
%   SHAPE = MEASURE_SPIKE_SHAPE_MEX(X,T,I) finds the peaks, troughs, and 
%   zero-crossings of the thresold-triggered spike waveforms specified by 
%   X,T,I.
%   
%   X is a MxCxN array of single-precision floating-point values. C is the 
%   number of recording channels in the C-trode, N is the number of 
%   trigger events, and M is the number of samples acquired per window.
%   
%   T is a Mx1 vector of single-precision floating-point values, which 
%   gives the times of the M samples relative to the first sample in the 
%   window. T is assumed to be monotonically increasing. 
%   
%   I is an integer scalar between 1 and M, which indexes along T and along 
%   the first dimension of X. I specifies the time of the trigger within 
%   the M-sample acquisition window. I is expressed with MATLAB-style 
%   (i.e., the first element of T corresponds to I = 1).
%   
%   The output argument SHAPE is a struct whose fields are NxC arrays of 
%   single-precision floating-point values"
%
%   Troughs and peaks are found by parabolic interpolation, and
%   rising/falling/crossing points are found by linear interpolation. The
%   interpolation can be improved by pre-processing with cubic spline
%   interpolation, which is what is done in the m-file MEASURE_SPIKE_SHAPE.
%
%Written by SMK, 2009 August 26.
