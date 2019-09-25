function n = length(a)
%LENGTH Length of a categorical array.
%   N = LENGTH(A), when A is not empty, returns the size of the longest
%   dimension of the categorical array A.  If A is a vector, this is the same
%   as its length.  LENGTH is equivalent to MAX(SIZE(X)) for non-empty arrays,
%   and 0 for empty arrays.
%
%   See also CATEGORICAL/NUMEL, CATEGORICAL/SIZE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:53 $

n = length(a.codes);
