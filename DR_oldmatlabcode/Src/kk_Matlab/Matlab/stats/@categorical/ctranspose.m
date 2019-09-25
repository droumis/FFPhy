function b = ctranspose(a)
%CTRANSPOSE Transpose a categorical matrix.
%   B = CTRANSPOSE(A) returns the transpose of the 2-dimensional categorical
%   matrix A.  Note that CTRANSPOSE is identical to TRANSPOSE for categorical
%   arrays.
%
%   CTRANSPOSE is called for the syntax A'.
%
%   See also CATEGORICAL/TRANSPOSE, CATEGORICAL/PERMUTE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:30:28 $

b = a;
b.codes = a.codes';
