function b = double(a,vars)
%DOUBLE Convert dataset variables to a double array.
%   B = DOUBLE(A) returns the contents of the dataset A, converted to one
%   double array.  The classes of the variables in the dataset must support
%   the conversion.
%
%   B = DOUBLE(A,VARS) returns the contents of the dataset variables specified
%   by VARS.  VARS is a positive integer, a vector of positive integers,
%   a variable name, a cell array containing one or more variable names, or a
%   logical vector.
%
%   See also DATASET, DATASET/SINGLE, DATASET/REPLACEDATA.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.8.2 $  $Date: 2006/11/11 22:56:35 $

if nargin < 2 || isempty(vars)
    vars = 1:a.nvars;
else
    vars = getvarindices(a,vars,false);
end

dims = cellfun('ndims',a.data(vars));
if any(diff(dims))
    error('stats:dataset:double:DimensionMismatch', ...
          'All dataset variables must have the same number of dimensions.');
end
sizes = cellfun(@size,a.data(vars),'uniformOutput',false);
sizes = cell2mat(sizes(:));
if any(any(diff(sizes(:,[1 3:end]),1),1))
    error('stats:dataset:double:SizeMismatch', ...
          'All dataset variables must have the same length in all but the second dimension.');
end

endCol = cumsum(sizes(:,2),1);
startCol = [1; endCol(1:end-1)+1];
szOut = sizes(1,:); szOut(2) = sum(sizes(:,2),1);
b = zeros(szOut);
for j = 1:length(vars)
    b(:,startCol(j):endCol(j),:) = double(a.data{vars(j)});
end
