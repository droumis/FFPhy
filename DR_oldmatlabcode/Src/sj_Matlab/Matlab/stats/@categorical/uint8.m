function b = uint8(a)
%UINT8 Convert categorical array to an uint8 array.
%   B = UINT8(A) converts the categorical array A to a uint8 array.  Each
%   element of B contains the internal categorical level code for the
%   corresponding element of A.
%
%   Undefined elements of A are assigned the value 0 in B.  If A contains
%   more than INTMAX('uint8') levels, the internal codes will saturate to
%   INTMAX('uint8') when cast to int8.
%
%   See also CATEGORICAL/DOUBLE, CATEGORICAL/INT8.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:20 $

b = uint8(a.codes);
