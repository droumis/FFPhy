function ai = interp1_angle(x,a,xi,jump_tol,varargin)
%INTERP1_ANGLE 1-D interpolation of angle.
%
%   AI = INTERP1_ANGLE(X,A,XI,JUMP_TOL) interpolates to find AI, the values of
%   the underlying function A, at the points in the array XI. X must be a
%   vector. JUMP_TOL specifies the tolerance for detecting wraparound
%   discontinuities in A, as described in the documentation for the MATLAB
%   UNWRAP function.
%
%   Angle values in A are assumed to be measured in radians on the interval
%   [-pi,+pi]. Any input values that are outside of this interval are
%   automatically wrapped into this interval.
%
%   If A is a vector, then it must also have length N, and the function will
%   return AI with the same size as XI. If A is an array, then it must be of
%   size [N,D1,D2,...Dk], and the interpolation is performed for each
%   D1-by-D2-by-...-Dk block of values in A(i,:,:,...,:).
%
%   If XI is a vector of length M, then AI has size [M,D1,D2,...,Dk]. If XI is
%   an array of size [M1,M2,...,Mj], then AI is of size
%   [M1,M2,...,Mj,D1,D2,...,Dk]
%
%   AI = INTERP1_ANGLE(X,A,XI,JUMP_TOL,METHOD) specifies alternate
%   interpolation methods. Documentation for the MATLAB function INTERP1
%   describes alternative methods. The default method of interpolation is
%   linear interpolation.
%
%   AI = INTERP1_ANGLE(X,A,XI,JUMP_TOL,METHOD,EXTRAP) specifies how to
%   extrapolate outside of the range of X as described in the documentation for
%   INTERP1.
%
%   See also INTERP1, UNWRAP.
%
%Depends on:
%   RECENTER_ANGLE (written by SMK)
%
%Written by smk 03 November 2008.
%

if (exist('recenter_angle') ~= 2)
  error('INTERP1_ANGLE depends on m-file RECENTER_ANGLE (written by SMK)');
end

% input checking
if ~isvector(x) || isempty(x) || ~isfloat(x) || ~isreal(x) || ~all(isfinite(x))
  error('X must be a non-empty vector of real floating-point finite values');
end
if isempty(a) || ~isfloat(a) || ~isreal(a)
  error('A must be a non-empty floating-point real array');
end
if (size(a,1) ~= numel(x))
  error('first dimension of A must be the same length as X');
end
if ~isequal(x,unique(x))
  error('values in X vector must be unique');
end
if any((a(:) <-pi) | (a(:) > pi))
  warning('values of A outside the interval [-pi,+pi] will be wrapped');
end
if any(isnan(a(:))) || any(~isfinite(a(:)))
  warning('interpolation near NaN or Inf elements of A will be undefined');
end
if ~isnumeric(xi) || ~isfloat(xi) || ~isreal(xi) || ~all(isfinite(xi(:)))
  error('XI must be a real floating-point array of finite values');
end
if ~isnumeric(jump_tol) || ~isscalar(jump_tol) || ~isreal(jump_tol) || ...
    ~isfinite(jump_tol) || ~isfloat(jump_tol)
  error('JUMP_TOL must be a real finite floating-point scalar');
end

% It's very important to perform all computations in double precision, because
% MATLAB's mod function accumulates floating-point errors
% Map values to the range [-pi,+pi] and then unwrap
a_unwrapped = unwrap(recenter_angle(double(a),0),jump_tol,1);
% call INTERP1 with varargin on the unwrapped angles and backtransform to
% recover values in [-pi,+pi]
ai = recenter_angle(interp1(x,a_unwrapped,xi,varargin{:}),0);

ai = cast(ai,class(a));

