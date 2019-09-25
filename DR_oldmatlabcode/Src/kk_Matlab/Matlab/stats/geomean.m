function m = geomean(x,dim)
%GEOMEAN Geometric mean.
%   M = GEOMEAN(X) returns the geometric mean of the values in X.  When X
%   is an n element vector, M is the n-th root of the product of the n
%   elements in X.  For a matrix input, M is a row vector containing the
%   geometric mean of each column of X.  For N-D arrays, GEOMEAN operates
%   along the first non-singleton dimension.
%
%   GEOMEAN(X,DIM) takes the geometric mean along dimension DIM of X.
%
%   See also MEAN, HARMMEAN, TRIMMEAN.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.10.2.3 $  $Date: 2004/07/28 04:38:23 $

if any(x(:) < 0)
    error('stats:geomean:BadData', 'X may not contain negative values.')
end

if nargin < 2 || isempty(dim)
    % Figure out which dimension sum will work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

n = size(x,dim);
% Prevent divideByZero warnings for empties, but still return a NaN result.
if n == 0, n = NaN; end

% Take the n-th root of the product of elements of X, along dimension DIM.
% prod(x).^(1/n) would not give logOfZero warnings, they're an artifact of
% using logs, so silence them.
warn = warning('off', 'MATLAB:log:logOfZero');
if nargin < 2
    m = exp(sum(log(x))./n);
else
    m = exp(sum(log(x),dim)./n);
end
warning(warn)
