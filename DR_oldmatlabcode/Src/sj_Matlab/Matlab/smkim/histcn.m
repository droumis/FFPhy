function [n, bin] = histcn(x,varargin);
%HISTCN Multi-dimensional histogram count
%
%   N = HISTCN(X,EDGES1,EDGES2,...,EDGESN), for matrix X, counts the number of
%   P-dimensional vector values (rows) in X that fall within the voxels
%   delimited by the vectors EDGES1,EDGES2,...,EDGESN. N is an array of size
%   CELLFUN(@LENGTH,{EDGES1, ..., EDGESN}). 
%
%   N(k1,k2,...,kN) will count the row X(i,:) if 
%       EDGES1(k1) <= X(i,1) < EDGES1(k1+1) &
%       EDGES2(k2) <= X(i,2) < EDGES2(k2+1) & ...
%       EDGESP(kP) <= X(i,P) < EDGESN(kP+1)
%   The last bins along each dimension of N will count any values that equal
%   EDGES*(end). Values outside of the range of values in EDGES* are not
%   counted. Use -Inf and +Inf to include all non-NaN values.
%
%   [N, BIN] = HISTCN(...) returns a matrix BIN of the same size as X, whose
%   rows are vectors of indices into the corresponding bins of EDGES1, EDGES2,
%   etc., where a value of 1 means that the value falls between the 1st and 2nd
%   edge value, a value of 2 means that the value falls between the 2nd and 3rd
%   edge value, ..., and a value of length(EDGES) means that the value equals
%   EDGES(end). A value of 0 means that the value did not fall within the range
%   of EDGES.
%
%Written by SMK, 2009 December 23.
%

if ~isfloat(x) || ~isreal(x)
  error('X must be floating-point real');
end 
if ndims(x) ~= 2
  error('X must be an M-by-N (2-dimensional) matrix');
end
% m is the number of observations
% p is the dimensionality of the observations
[m, p] = size(x);
% edges{1} defines histogram bins along the 1st dimension, 
% edges{2} defines bins along the 2nd dimension, etc.
edges = varargin;
% Check that dimensionality of observations equals the 
% number of EDGES* vectors provided
if (numel(edges) ~= p)
  error('Number of EDGES vectors must equal number of columns in X');
end
if ~all(cellfun(@(e) isreal(e) && isvector(e) && all(diff(e) > 0),edges))
  error('EDGES vectors must be real numeric and monotonic increasing');
end

% Allocate output array which has the number of bins 
% jointly defined by the EDGES* vectors
sz = cellfun(@length,edges);
% For degenerate case of one-dimensional histogram, ensure that output is a
% column vector
if (p == 1)
  sz = [sz 1];
end
% For each vector observation X(i,:), BIN(i,:) summarizes where that vector
% falls into the bins defined by EDGES.
bin = zeros(m,p);
for i = 1:p
  % Use HISTC to count observations in bins defined by EDGES{i}. Note that,
  % because of how HISTC works, observations which do not fall within the
  % specified edges are assigned zero values in BIN(:,i)
  [junk, bin(:,i)] = histc(x(:,i),edges{i});
end
% Accumulate counts for rows of BIN that in which all subscripts are nonzero
nzsubs = find(all(bin,2));
n = accumarray(bin(nzsubs,:),ones([size(nzsubs,1) 1]),sz,@sum,0);

