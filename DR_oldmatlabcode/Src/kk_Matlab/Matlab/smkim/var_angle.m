function v = var_angle(a,w,dim)
%VAR_ANGLE Circular variance.
%   V = VAR_ANGLE(A) returns the circular variance of the angle values in A.
%   Values in A are in radians. If A is a matrix, V is a row vector containing
%   the variance of each column in A. For N-D arrays, VAR_ANGLE operates along
%   the first non-singleton dimension of A.
%
%   V = VAR_ANGLE(A,W) computes the circular variance using the weight vector W.
%   The length of W must equal the length of the dimension over which VAR_ANGLE
%   operates, and its elements must be nonnegative. VAR_ANGLE normalizes W to
%   sum to one. If W is [], no special weighting is applied.
%
%   V = VAR_ANGLE(A,W,DIM) takes the circular variance along the dimension DIM
%   of X.  Pass in [] for W to weight all elements of A equally.
%
%   See also VAR.
%
%Written by smk, 3 November 2008.
%

if ~isfloat(a) || ~isreal(a) || ~all(isfinite(a(:)))
  error('A must be floating-point real with no Inf values');
end
if any((a(:) < -pi) | (a(:) > pi))
  warning('some A values are outside the interval [-pi,+pi]');
end
if nargin<2
  w = [];
end
if (nargin < 3)
  % If A is [], return NaN
  if isequal(a,[])
    v = nan(class(a)); 
    return; 
  end
  % Identify first non-singleton dimension of A
  dim = find(size(a)~=1,1);
  if isempty(dim)
    dim = 1; 
  end
end
n = size(a,dim);

if isempty(w)
  z = complex(cos(a),sin(a));
else
  if ~isvector(w) || any(w < 0)
    error('W must be a vector of non-negative weights');
  end
  if (numel(w) ~= n)
    error(['length of W must match the size of A along the specified ' ...
        'dimension DIM']);
  end    
  % normalize W, and reshape and tile to align with A
  w = w ./ mean(w);
  wresize = ones([1 max(ndims(a),dim)]); 
  wresize(dim) = n;
  wtile = size(a); 
  if dim <= ndims(a)
    wtile(dim) = 1; 
  end
  w = repmat(reshape(w,wresize),wtile);
  % element-wise weighting of values in A
  assert(abs(mean(w) - 1) < sqrt(eps(max(w))));
  z = w .* complex(cos(a),sin(a));
end
v = 1 - abs(mean(z,dim));



