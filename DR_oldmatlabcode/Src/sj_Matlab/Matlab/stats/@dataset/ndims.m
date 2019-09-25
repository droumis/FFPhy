function n = ndims(a)
%NDIMS Number of dimensions of a dataset array.
%   N = NDIMS(A) returns the number of dimensions in the dataset A.  The
%   number of dimensions in an array is always 2.
%  
%   See also DATASET/SIZE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:43 $

n = a.ndims;
