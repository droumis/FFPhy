function b = int32(a)
%INT32 Convert categorical array to an int32 array.
%   B = INT32(A) converts the categorical array A to an int32 array.  Each
%   element of B contains the internal categorical level code for the
%   corresponding element of A.
%
%   Undefined elements of A are assigned the value 0 in B.
%
%   See also CATEGORICAL/DOUBLE, CATEGORICAL/UINT32.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:41 $

b = int32(a.codes);
