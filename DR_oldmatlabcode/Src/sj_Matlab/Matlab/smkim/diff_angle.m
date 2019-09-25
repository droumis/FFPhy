function d = diff_angle(a,jump_tol,varargin);
%DIFF_ANGLE Difference and approximate derivatives for angles.
%
%   DIFF_ANGLE works like DIFF, except that it correctly handles the wraparound
%   of circular data. 
%
%   DIFF_ANGLE(A,JUMP_TOL), for a vector A of angles (in radians), is the
%   vector of consecutive differences. JUMP_TOL specifies the tolerance for
%   detecting wraparound discontinuities in A, as described in the
%   documentation for the MATLAB UNWRAP function.
%
%   DIFF_ANGLE(A,JUMP_TOL), for a matrix A, is the matrix of row differences
%   (i.e., the same result as appplying DIFF_ANGLE separately to each column of
%   A and horizontally concatenating the results).
%
%   DIFF_ANGLE(A,JUMP_TOL), for an N-D array A, is the difference along the
%   first non-singleton dimension of A.
%
%   DIFF_ANGLE(A,JUMP_TOL,N) is the N-th order difference along the first
%   non-singleton dimension (denote it by DIM). If N >= size(X,DIM), DIFF_ANGLE
%   takes successive differences along the next non-singleton dimension.
%
%   DIFF_ANGLE(X,JUMP_TOL,N,DIM) is the Nth difference function along dimension
%   DIM. If N >= size(X,DIM), DIFF_ANGLE returns an empty array.
%
%   See also DIFF, UNWRAP.
%
%Depends on:
%   RECENTER_ANGLE (written by SMK)
%
%Written by smk 03 November 2008.
%

if (exist('recenter_angle') ~= 2)
  error('DIFF_ANGLE depends on m-file RECENTER_ANGLE (written by SMK)');
end

if ~isfloat(a) || ~isreal(a)
  error('A must be floating-point real');
end
if any((a < -pi) | (a > pi))
  %warning('some A values are outside the interval [-pi,+pi]');
end
if ~isscalar(jump_tol) || ~isfloat(jump_tol) || ~isreal(jump_tol) || ...
    ~isfinite(jump_tol)
  error('JUMP_TOL must be a real finite floating-point scalar');
end

% Unwrap along the specified dimension.
if (length(varargin) == 2)
  dim = varargin{2};
else
  dim = 1;
end
% It's very important to perform all computations in double precision, because
% MATLAB's mod function accumulates floating-point errors
a_unwrapped = unwrap(recenter_angle(double(a),0),jump_tol,dim);

% Compute diff and remap to standard interval [-pi,+pi]
d = recenter_angle(diff(a_unwrapped,varargin{:}),0);

d = cast(d,class(a));

