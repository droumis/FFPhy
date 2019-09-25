function b = single(a)
%SINGLE Convert categorical array to single array.
%   B = SINGLE(A) converts the categorical array A to a single array.  Each
%   element of B contains the internal categorical level code for the
%   corresponding element of A.
%
%   Undefined elements of A are assigned the value NaN in B.
%
%   See also CATEGORICAL/DOUBLE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:09 $

b = single(a.codes);
b(b==0) = NaN;