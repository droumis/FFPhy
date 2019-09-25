function b = flipud(a)
%FLIPUD Flip categorical matrix in up/down direction.
%   B = FLIPUD(A) returns the 2-dimensional categorical matrix A with columns
%   preserved and rows flipped in the up/down direction.
%
%   See also CATEGORICAL/FLIPDIM,  CATEGORICAL/FLIPUD, CATEGORICAL/ROT90.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:37 $

b = a;
b.codes = flipud(a.codes);
