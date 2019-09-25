function t = isempty(a)
%ISEMPTY True for empty dataset array.
%   TF = ISEMPTY(A) returns true (1) if A is an empty dataset and false (0)
%   otherwise. An empty array has no elements, that is PROD(SIZE(A))==0.
%  
%   See also DATASET/SIZE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/12/15 19:31:40 $

t = (a.nobs == 0) || (a.nvars == 0);
