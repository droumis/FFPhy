function zi = interp2_torus(a,b,z,ai,bi,varargin)
%INTERP2_TORUS 2-D interpolation on a torus.
%
%   ZI = INTERP2_TORUS(A,B,Z,AI,BI) interpolates to find ZI, the values of the
%   underlying 2-D function Z at the points in matrices AI and BI. Matrices
%   A and B specify the points on a toroidal manifold at which the data Z is
%   given. A and B correspond to the periodic dimensions of the torus. All
%   elements in A and B must fall within the interval [-pi,+pi].
%
%   AI can be a row vector, in which case it specifies a matrix with constant
%   columns. Similarly, BI can be a column vector and it specifies a matrix with
%   constant rows.
%
%   ZI = INTERP2_CYLINDER(...,METHOD) specifies alternate methods. The default
%   is linear interpolation. Available methods are described in the dcumentation
%   for the MATLAB function INTERP2.
%
%   All the interpolation methods require that X and A be monotonic and
%   plaid (as if they were created using MESHGRID).  If you provide two
%   monotonic vectors, INTERP2_CYLINDER changes them to a plaid internally.
%   X and A can be non-uniformly spaced.
%
%Written by SMK, 2010 January 28.
%

% Check correctness of A
if isempty(a) || (ndims(a) > 2) || ~isfloat(a) || ~isreal(a) || ...
    ~all(isfinite(a(:))) || any((a(:) < -pi) | (a(:) > +pi))
  error(['A must be a non-empty vector or matrix containing ' ...
      'floating-point real values on the interval [-pi,+pi]']);
end
% If A is a vector, it must be a row vector
if isvector(a) && (size(a,1) > size(a,2))
  error('If A is a vector, it must be a row vector');
end
tf = bsxfun(@eq,a(1,:),a);
if ~all(tf(:))
  error(['A must be a vector, or a plaid matrix with identical ' ...
    'rows of the type produced by MESHGRID(A,B)']);
end
sn = unique(sign(diff(a)));
if (numel(sn) > 1) || (sn == 0) || ...
    ((sn > 0) && any(a(:,1) + 2*pi <= a(:,end))) || ...
    ((sn < 0) && any(a(:,1) - 2*pi >= a(:,end)))
  error(['Rows of A must be monotonic after unwrapped modulo (2*pi)']);
end

% Check correctness of B
if isempty(b) || (ndims(b) > 2) || ~isfloat(b) || ~isreal(b) || ...
    ~all(isfinite(b(:))) || any((b(:) < -pi) | (b(:) > +pi))
  error(['B must be a non-empty vector or matrix containing ' ...
      'floating-point real values on the interval [-pi,+pi]']);
end
% If B is a vector, it must be a column vector
if isvector(b) && (size(b,1) < size(b,2))
  error('If B is a vector, it must be a column vector');
end
tf = bsxfun(@eq,b(:,1),b);
if ~all(tf(:))
  error(['B must be a vector, or a plaid matrix with identical ' ...
    'columns of the type produced by MESHGRID(A,B)']);
end
sn = unique(sign(diff(b)));
if (numel(sn) > 1) || (sn == 0) || ...
    ((sn > 0) && any(b(1,:) + 2*pi <= b(end,:))) || ...
    ((sn < 0) && any(b(1,:) - 2*pi >= b(end,:)))
  error(['Columns of B must be monotonic after unwrapped modulo (2*pi)']);
end

% AI and BI are allowed to contain NaN values
if ~isfloat(ai) || ~isreal(ai) || any(isinf(ai(:)))
  error(['AI must contain finite floating-point real values']);
end
if any((ai(:) < -pi) | (ai(:) > +pi))
  warning(['Some elements of AI lie utside the interval [-pi,+pi]. ' ...
      'These will be automatically wrapped']);
end
if ~isfloat(bi) || ~isreal(bi) || any(isinf(bi(:)))
  error(['BI must contain finite floating-point real values']);
end
if any((bi(:) < -pi) | (bi(:) > +pi))
  warning(['Some elements of BI lie utside the interval [-pi,+pi]. ' ...
      'These will be automatically wrapped']);
end

% Tile the data to cover periodic boundaries.
a = repmat([a-2*pi, a, a+2*pi],[3 1]);
b = repmat([b-2*pi; b; b+2*pi],[1 3]);
z = repmat(z,[3 3]);

try
  zi = interp2(a,b,z,ai,bi,varargin{:});
catch
  error('Error while calling INTERP2. Perhaps input arguments are not valid?');
end

