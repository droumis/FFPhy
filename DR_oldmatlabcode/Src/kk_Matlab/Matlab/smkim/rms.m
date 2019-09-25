function r = rms(x,m)
%RMS Compute root-mean-square of a real vector
%
%   R = RMS(X,M) takes a floating-point real vector X containing finite non-NaN
%   values and a floating-point real scalar M, and returns a scalar value R
%   which is the root-mean-square of X with respect to baseline M.
%
%   RMS(X) is the same as RMS(X,0).
%
%Written by SMK, 2009 November 18.
%

if (nargin == 1)
  m = 0;
end

if ~isvector(x) || ~isfloat(x) || isempty(x) || ~isreal(x) || ...
    ~all(isfinite(x))
  error('X must be a non-empty real floating-point vector');
end

if ~isscalar(m) || ~isfloat(m) || ~isreal(m) || ~isfinite(m)
  error('M must be a finite real floating-point scalar');
end

n = numel(x);
r = sqrt(sum((x - m).^2) / n);

