function dim = isvector(vect)
%ISVECTOR          True for 1-D arrays.
%   ISVECTOR(VECT) is non-zero if exactly one dimension of VECT has length
%   greater than 1.  The return value is then the index of that dimension.
%   Note that NDIMS can not be used to decide this question, because it
%   returns 2 for, e.g., (M x 1) and (1 x M) arrays.
%
%   Example:
%      isvector(1);             % returns 0
%      isvector([1 2 ; 3 4])    % returns 0
%      isvector([1:10])         % returns 2

nonsingle = [size(vect) > 1];
dim = find(nonsingle);
if ((length(dim)>1) || (isempty(dim))),  dim = 0;  end;