function s = std_angle(varargin)
%STD_ANGLE Circular standard deviation.
%   S = STD_ANGLE(A) returns the circular standard deviation of the angle values
%   in A. Note that the circular standard deviation is not just the square root
%   of the circular variance! 
%
%   Values in A are in radians. If A is a matrix, S is a row vector containing
%   the circular standard deviation of each column in A. For N-D arrays,
%   STD_ANGLE operates along the first non-singleton dimension of A.
%
%   S = STD_ANGLE(A,W) computes the circular standard deviation using the weight
%   vector W. The length of W must equal the length of the dimension over which
%   STD_ANGLE operates, and its elements must be nonnegative. STD_ANGLE
%   normalizes W to sum to one. If W is [], no special weighting is applied.
%
%   S = STD_ANGLE(A,W,DIM) takes the circular variance along the dimension DIM
%   of X.  Pass in [] for W to weight all elements of A equally.
%
%   See also STD.
%
%Depends on:
%   VAR_ANGLE (written by smk)
%
%Written by smk, 3 November 2008.
%

if (exist('var_angle') ~= 2)
  error('STD_ANGLE depends on m-file VAR_ANGLE (written by smk)');
end

v = var_angle(varargin{:});
s = sqrt(-2 * log(1-v));

