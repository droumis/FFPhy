function b = uint16(a)
%UINT16 Convert categorical array to an uint16 array.
%   B = UINT16(A) converts the categorical array A to a uint16 array.  Each
%   element of B contains the internal categorical level code for the
%   corresponding element of A.
%
%   Undefined elements of A are assigned the value 0 in B.
%
%   See also CATEGORICAL/DOUBLE, CATEGORICAL/INT16.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:17 $

b = uint16(a.codes);
