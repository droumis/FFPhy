function m = mean_angle(a,dim)
%MEAN_ANGLE Mean direction for circular (angle) data.
%   M = MEAN_ANGLE(A) returns the circular mean of the angle values in A. Values
%   in A are in radians. If A is a matrix, M is a row vector containing the
%   means of the columns of A. For N-D arrays, MEAN_ANGLE operates along the
%   first non-singleton dimension of A.
%
%   V = MEAN_ANGLE(A,DIM) takes the mean along dimension DIM.
%
%   See also MEAN.
%
% Written by smk, 3 November 2008.
%

if ~isfloat(a) || ~isreal(a) || ~all(isfinite(a(:)))
  error('A must be floating-point real with no Inf values');
end
if any((a < -pi) | (a > pi))
  warning('some A values are outside the interval [-pi,+pi]');
end

if (nargin == 1)
  % Determine which dimension to take mean along
  dim = min(find(size(a)~=1));
  if isempty(dim)
    dim = 1;
  end
end
z = sum(complex(cos(a),sin(a)),dim);
m = angle(z);
% Circular mean is undefined if resultant vector is zero
m(z == complex(0,0)) = NaN;

