function b = uint32(a)
%UINT32 Convert categorical array to an uint32 array.
%   B = UINT32(A) converts the categorical array A to a uint32 array.  Each
%   element of B contains the internal categorical level code for the
%   corresponding element of A.
%
%   Undefined elements of A are assigned the value 0 in B.
%
%   See also CATEGORICAL/DOUBLE, CATEGORICAL/INT32.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:18 $

b = uint32(a.codes);
