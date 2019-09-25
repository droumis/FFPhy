function b = double(a)
%DOUBLE Convert categorical array to double array.
%   B = DOUBLE(A) converts the categorical array A to a double array.  Each
%   element of B contains the internal categorical level code for the
%   corresponding element of A.
%
%   Undefined elements of A are assigned the value NaN in B.
%
%   See also CATEGORICAL/SINGLE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:31 $

b = double(a.codes);
b(b==0) = NaN;