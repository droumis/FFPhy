function e = end(a,k,n)
%END Last index in an indexing expression for a dataset array.
%   END(A,K,N) is called for indexing expressions involving the dataset A
%   when END is part of the K-th index out of N indices.  For example, the
%   expression A(end-1,:) calls A's END method with END(A,1,2).
%
%   See also DATASET/SIZE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.3 $  $Date: 2006/12/15 19:31:36 $

switch k
case 1
    e = a.nobs;
case 2
    e = a.nvars;
otherwise
    e = 1;
end
