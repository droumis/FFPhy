function b = int8(a)
%INT8 Convert categorical array to an int8 array.
%   B = INT8(A) converts the categorical array A to an int8 array.  Each
%   element of B contains the internal categorical level code for the
%   corresponding element of A.
%
%   Undefined elements of A are assigned the value 0 in B.  If A contains
%   more than INTMAX('int8') levels, the internal codes will saturate to
%   INTMAX('int8') when cast to int8.
%
%   See also CATEGORICAL/DOUBLE, CATEGORICAL/UINT8.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:43 $

b = int8(a.codes);
