function b = int16(a)
%INT16 Convert categorical array to an int16 array.
%   B = INT16(A) converts the categorical array A to an int16 array.  Each
%   element of B contains the internal categorical level code for the
%   corresponding element of A.
%
%   Undefined elements of A are assigned the value 0 in B.
%
%   See also CATEGORICAL/DOUBLE, CATEGORICAL/UINT16.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:40 $

b = int16(a.codes);
