function zi = interp2_cylinder(x,a,z,xi,ai,varargin)
%INTERP2_CYLINDER 2-D interpolation on a cylinder.
%
%   ZI = INTERP2_CYLINDER(X,A,Z,XI,AI) interpolates to find ZI, the values of
%   the underlying 2-D function Z at the points in matrices XI and AI. Matrices
%   X and A specify the points on a cylindrical manifold at which the data Z
%   is given. X corresponds to the longitudinal axis of the cylinder, and A
%   corresponds to the wrapped circumference of the cylinder. All elements in A
%   must fall within the interval [-pi,+pi].
%
%   XI can be a row vector, in which case it specifies a matrix with constant
%   columns. Similarly, AI can be a column vector and it specifies a matrix with
%   constant rows.
%
%   ZI = INTERP2_CYLINDER(...,METHOD) specifies alternate methods. The default
%   is linear interpolation. Available methods are described in the dcumentation
%   for the MATLAB function INTERP2.
%
%   ZI = INTERP2_CYLINDER(...,METHOD,EXTRAPVAL) specificies a method and a
%   scalar value for ZI outside of the range of X. Thus, ZI will equal EXTRAPVAL
%   for any value of XI which is not spanned by X. A method must be specified
%   for EXTRAPVAL to be used, the default method is 'linear'.
%
%   All the interpolation methods require that X and A be monotonic and
%   plaid (as if they were created using MESHGRID).  If you provide two
%   monotonic vectors, INTERP2_CYLINDER changes them to a plaid internally.
%   X and A can be non-uniformly spaced.
%
%Written by SMK, 2010 January 28.
%


% (We allow INTERP2 take care of checking X)

% Check correctness of A
if isempty(a) || (ndims(a) > 2) || ~isfloat(a) || ~isreal(a) || ...
    ~all(isfinite(a(:))) || any((a(:) < -pi) | (a(:) > +pi))
  error(['A must be a non-empty vector or matrix containing ' ...
      'floating-point real values on the interval [-pi,+pi]']);
end
% If A is a vector, it must be a column vector
if isvector(a) && (size(a,1) < size(a,2))
  error('If A is a vector, it must be a column vector');
end
tf = bsxfun(@eq,a(:,1),a);
if ~all(tf(:))
  error(['A must be a vector, or a plaid matrix with identical ' ...
    'columns of the type produced by MESHGRID(X,A)']);
end
sn = unique(sign(diff(a)));
if (numel(sn) > 1) || (sn == 0) || ...
    ((sn > 0) && any(a(1,:) + 2*pi <= a(end,:))) || ...
    ((sn < 0) && any(a(1,:) - 2*pi >= a(end,:)))
  error(['Columns of A must be monotonic after unwrapped modulo (2*pi)']);
end

% AI is allowed to contain NaN values
if ~isfloat(ai) || ~isreal(ai) || any(isinf(ai(:)))
  error(['AI must contain finite floating-point real values']);
end
if any((ai(:) < -pi) | (ai(:) > +pi))
  warning(['Some elements of AI lie utside the interval [-pi,+pi]. ' ...
      'These will be automatically wrapped']);
end

% Tile the data to cover periodic boundaries
a = [a-2*pi; a; a+2*pi];
z = repmat(z,[3 1]);
if ~isvector(x)
  % If X is a matrix, then we also need to tile X
  x = repmat(x,[3 1]);
end

try
  zi = interp2(x,a,z,xi,ai,varargin{:});
catch
  error('Error while calling INTERP2. Perhaps input arguments are not valid?');
end


