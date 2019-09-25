function b = uint64(a)
%UINT64 Convert categorical array to an uint64 array.
%   B = UINT64(A) converts the categorical array A to a uint64 array.  Each
%   element of B contains the internal categorical level code for the
%   corresponding element of A.
%
%   Undefined elements of A are assigned the value 0 in B.
%
%   See also CATEGORICAL/DOUBLE, CATEGORICAL/INT64.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:19 $

b = uint64(a.codes);
