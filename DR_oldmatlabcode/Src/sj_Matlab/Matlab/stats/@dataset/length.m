function len = length(a)
%LENGTH Length of a dataset array.
%   N = LENGTH(A) returns the number of observations in the dataset A.  LENGTH
%   is equivalent to SIZE(A,1).
%  
%   See also DATASET/SIZE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:42 $

len = a.nobs;
