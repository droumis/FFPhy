function t = isscalar(a)
%ISSCALAR True if categorical array is a scalar.
%   TF = ISSCALAR(A) returns true (1) if the categorical array A is a 1-by-1
%   matrix, and false (0) otherwise.
%
%   See also CATEGORICAL/ISVECTOR, CATEGORICAL/ISEMPTY, CATEGORICAL/SIZE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:50 $

t = isscalar(a.codes);
