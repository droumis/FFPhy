function b = squeeze(a)
%SQUEEZE Squeeze singleton dimensions from a categorical array.
%   B = SQUEEZE(A) returns an array B with the same elements as the
%   categorical array A but with all the singleton dimensions removed.  A
%   singleton is a dimension such that size(A,DIM)==1.  2-D arrays are
%   unaffected by SQUEEZE so that row vectors remain rows.
%
%   See also CATEGORICAL/SHIFTDIM.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:11 $

b = a;
b.codes = squeeze(a.codes);
