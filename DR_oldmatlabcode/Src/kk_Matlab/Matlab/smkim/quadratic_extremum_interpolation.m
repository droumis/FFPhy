function [xi, yi] = quadratic_extremum_interpolation(x,y)
%QUADRATIC_EXTREMUM_INTERPOLATION Interpolate 1-dimensional extrema given neighboring three points
%
%   [XI, YI] = QUADRATIC_EXTREMUM_INTERPOLATION(X,Y) approximates local minima
%   and maxima of a discretely-sampled smooth univariate real scalar-valued
%   function. X and Y are N-by-3 matrices of real floating-point values, where
%   each row specifies a triplet of (x,y) points that are sampled around a local
%   minimum or maximum of a function f: x -> y. The function f is assumed to
%   have continuous second derivative with respect to x, so that it can be
%   locally approximated by a quadratic function. The inputs must meet the
%   following conditions:
%
%     TOL = max(eps(X),[],2)
%     all((X(:,2) - X(:,1) > TOL) & (X(:,3) - X(:,2) > TOL))
%     all(sign(Y(:,2) - Y(:,1)) ~= sign(Y(:,3) - Y(:,2)))
%
%
%   XI and YI are N-by-1 column vectors of type double.
%
%Written by SMK, 2010 January 5.
%

if ~isfloat(x) || ~isfloat(y) || ~isreal(x) || ~isreal(y) || ...
    ~isequal(size(x),size(y)) || ~isequal(2,ndims(x),ndims(x)) || ...
    (size(x,2) ~= 3) || ...
    ~all(all(bsxfun(@gt,diff(x,1,2),max(eps(x),[],2)),2),1) || ...
    ~all(sign(y(:,2) - y(:,1)) ~= sign(y(:,3) - y(:,2)))
  error('Inputs are not valid');
end

% Recentered values in each row for numerical stability
dx12 = x(:,1) - x(:,2);
dx32 = x(:,3) - x(:,2);
dx31 = x(:,3) - x(:,1);

dy12 = y(:,1) - y(:,2);
dy32 = y(:,3) - y(:,2);
dy31 = y(:,3) - y(:,1);

c = dx12 .* dx32 .* dx31;
a = (dx12.*dy32 - dx32.*dy12)./c;
b = (dy12.*(dx32.^2) - dy32.*(dx12.^2))./c;

xi = x(:,2) - b./(2*a);
yi = y(:,2) - (b.^2)./(4*a);

