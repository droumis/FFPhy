function y = find_blobs(x,periodic,thresholds)
%FIND_BLOBS Segment contiguous blobs in an N-dimensional gridded intensity map using watershed algorithm.
%
%   Y = FIND_BLOBS(X,PERIODIC,THRESHOLDS) finds contiguous "blobs" around
%   local maxima in the gridded data X using a watershed algorithm. 
%
%   X must be an N-dimensional array whose elements are all finite. Values in X
%   are assumed to reside on a uniformly-spaced grid. 
%
%   PERIODIC must be a logical row vector whose length equals the number of
%   dimensions of X. If PERIODIC(d) is true, then X is interpreted as being
%   periodic along its dth dimension. (For example, if X is a 2-dimensional
%   array, and PERIODIC = [true true], then X is a torus.)
%
%   THRESHOLDS is a two-element vector of threshold values of the same numeric
%   class as X. A blob is defined as a collection of contiguous voxels whose
%   values are all greater than or equal to THRESHOLDS(1), such that at least
%   one of them has a value greater than THRESHOLDS(2). A pair of voxels are
%   "contiguous" if they share a common edge (for 2-dimensional data, this is
%   commonly referred to as "4-connectivity" in the morphological image
%   processing literature). Voxels whose values are greater than or equal to
%   THRESHOLDS(1) are partitioned into watersheds using local gradients.
%
%   The output Y is an array of the same size as X, containing non-negative
%   integer values. Voxels in X which do not belong to a blob are indicated by
%   zero. Voxels which belong to a blob are indicated by an integer value which
%   enumerates blob identity. The order of enumeration is arbitrary.
%   
%Depends on:
%   IMIMPOSEMIN (MATLAB Image Processing Toolbox)
%
%Written by smk, 2008 November 8.
%

if (exist('imimposemin') ~= 2)
  error(['FIND_BLOBS depends on the m-file IMIMPOSEMIN ' ...
      '(MATLAB Image Processing Toolbox)']);
end

if ~isreal(thresholds) || ~isvector(thresholds) || ...
    (numel(thresholds) ~= 2) || ~all(isfinite(thresholds))
  error('THRESHOLDS must be a two-element real vector with finite values');
end
if ~(thresholds(1) <= thresholds(2))
  error('THRESHOLDS(1) must be less than or equal to THRESHOLDS(2)');
end

if isempty(x) || ~isnumeric(x) || ~isreal(x) || ~all(isfinite(x(:)))
  error(['X must be a non-empty real numeric array of finite values']);
end
% Get size of X
sz = size(x);
dim = numel(sz);
% Subscripts and linear indices of voxels
num_voxels = prod(sz);
ind = reshape(1:num_voxels,sz);
[subs{1:dim}] = ind2sub(sz,ind);

if ~islogical(periodic) || ~isvector(periodic) || ...
    ~isequal(size(periodic),size(sz))
  error(['PERIODIC must be a logical row vector of size [1, ndims(X)] ' ...
      '(remember that MATLAB gives column vectors a dimensionality of 2)']);
end

% Apply IMIMPOSEMIN to force contiguous voxels whose values are less than
% thresholds(1) to be -Inf. To correctly respect periodic contiguity, we need to
% tile X along its periodic dimensions, apply IMIMPOSEMIN, and then snip out the
% central portion corresponding to the original array.
x = repmat(x,2*periodic + 1);
for i = 1:dim
  middle_subs{i} = subs{i} + sz(i)*periodic(i);
end
% IMIMPOSEMIN with (2*dim) connectivity
x = imimposemin(x,x < thresholds(1),2*dim);
x = x(sub2ind(sz .* (1+periodic),middle_subs{:}));
candidates = (x > -Inf);

% Sort voxels according to their adjusted values
[sort_val, sort_idx] = sort(x(:),1,'descend');
reverse_sort_idx = ind(sort_idx);

% For each voxel, get the linear indices of its neighbors which have values
% greater than -Inf; for each voxel, the maximum possible number of qualifying
% neighbors is (2*dim)
neighbors = cell(sz);
for i = 1:dim
  % Skip singleton dimension of X
  if (sz(i) == 1)
    continue;
  end
  flag = (1:dim == i);
  % Circularly shift subscripts to compute adjacency
  prev = circshift(ind,-1*flag);
  next = circshift(ind,+1*flag);
  if ~periodic(i)
    prev(subs{i} == 1) = 0;
    next(subs{i} == sz(i)) = 0;
  end
  % Exclude -Inf voxels
  next(x(next) == -Inf) = 0;
  prev(x(prev) == -Inf) = 0;
  % Zero is a sentinel index to indicate a missing or invalid neighbor
  neighbors = cellfun(@(c1,c2) cat(1,c1,c2),neighbors, ...
      arrayfun(@(x1,x2) cat(1,x1(x1 ~= 0),x2(x2 ~= 0)),prev,next, ...
      'UniformOutput',false),'UniformOutput',false);
end



% Identify voxels which are above thresholds(2) and which do not have any
% lower-valued neighbors; these are used as starting seeds
seeds = find( (x >= thresholds(2)) & ...
    arrayfun(@(i) all(x(neighbors{i}) <= x(i)),ind) );

%{
Lindeberg's watershed-based gray-level blob detection algorithm.

For simplicity, let us consider the case of detecting bright grey-level blobs
and let the notation "higher neighbour" stand for "neighbour pixel having a
higher grey-level value". Then, at any stage in the algorithm (carried out in
decreasing order of intensity values) is based on the following classification
rules:

[1] If a region has no higher neighbour, then it is a local maximum and will be the
    seed of a blob.
[2] Else, if it has at least one higher neighbour, which is background, then it
    cannot be part of any blob and must be background.
[3] Else, if it has more than one higher neighbour and if those higher
    neighbours are parts of different blobs, then it cannot be a part of any
    blob, and must be background.
[4] Else, it has one or more higher neighbours, which are all parts of the same
    blob. Then, it must also be a part of that blob.
%}

% Initialize output
y = zeros(sz);
% Keep track of which voxels have been checked
checked = (x == -Inf);
% Integer label to enumerate distinct blobs
label = 0;

% Function which returns Inf as the "maximum" of an empty matrix
function y = conditional_max(x)
  if isempty(x)
    y = -Inf;
  else
    y = max(x);
  end

