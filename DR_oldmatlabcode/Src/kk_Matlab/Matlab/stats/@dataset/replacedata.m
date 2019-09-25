function b = replacedata(a,x,vars)
%REPLACEDATA Convert array to dataset variables.
%   B = REPLACEDATA(A,X) creates a dataset B with the same variables as the
%   dataset A, but with the data for those variables replaced by the data in
%   the array X.  REPLACEDATA creates each variable in B using one or more
%   columns from X, in order.  X must have as many columns as the total number
%   of columns in all the variables in A, and as many rows as A has
%   observations.
%
%   B = REPLACEDATA(A,X,VARS) creates a dataset B with the same variables as
%   the dataset A, but with the data for the variables specified in VARS
%   replaced by the data in the array X. The remaining variables in B are
%   simply copies of the corresponding variables in A.  VARS is a positive
%   integer, a vector of positive integers, a variable name, a cell array
%   containing one or more variable names, or a logical vector.  Each variable
%   in B has as many columns as the corresponding variable in A.  X must have
%   as many columns as the total number of columns in all the variables
%   specified in VARS.
%
%   See also DATASET, DATASET/DOUBLE, DATASET/SINGLE.

%   Copyright 2006 The MathWorks, Inc. 
%   $Revision: 1.1.6.1 $  $Date: 2006/11/11 22:56:41 $

if nargin < 3 || isempty(vars)
    vars = 1:a.nvars;
else
    vars = getvarindices(a,vars,false);
end

szX = size(x);
ncols = cellfun('size',a.data(vars),2);

if szX(1) ~= a.nobs
    error('stats:dataset:replacedata:DimensionMismatch', ...
          'X must have the same number of rows as A has observations.');
elseif szX(2) ~= sum(ncols)
    error('stats:dataset:replacedata:DimensionMismatch', ...
          'The variables being replaced in A must have the same total number of columns as X.');
end

szOut = szX;
endCol = cumsum(ncols);
startCol = [1 endCol(1:end-1)+1];
b = a;
for j = 1:length(vars)
    szOut(2) = endCol(j) - startCol(j) + 1;
    b.data{vars(j)} = reshape(x(:,startCol(j):endCol(j),:), szOut);
end
