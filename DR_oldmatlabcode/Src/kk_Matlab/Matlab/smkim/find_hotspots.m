function h = find_hotspots(d,dims,thresh)
%FIND_HOTSPOTS Find regions of high density in a N-D density map
%
%   H = FIND_HOTSPOTS(D,DIMS,THRESH) returns an array of the same size as D,
%   containing nonnegative integer values. Zero values indicate that the
%   corresponding voxel in D does not belong to a hotspot. Voxels which are part
%   of a hotspot are labeled with enum integer values 1, 2, 3, ...
%   
%   D is a kernel density estimate obtained by DENSITYMAP
%
%   DIMS is a struct array of length N, with the following fields:
%       NAME is a human-readable descriptive string (e.g.'distance (cm)').
%       TYPE is 'linear' or 'circular'
%       GRID is a vector of evenly-spaced grid points at which the kernel
%       density is to be estimated; if TYPE is 'circular', then GRID(1)
%       coincides with GRID(end)
%
%   THRESH is a two-element vector of threshold values, expressed in the same
%   units as D. A hotspot is defined as a collection of contiguous voxels, of
%   which at least one voxel exceeds TRESH(1) and all voxels exceed THRESH(2).
%
%Written by smk, 2008 November 8.
%

if ~isnumeric(thresh) || ~isreal(thresh) || ~isvector(thresh) || ...
    (numel(thresh) ~= 2) || any(thresh > 1) || any(thresh < 0)
  error('THRESH must be a two-element vector with values between zero and 1');
end
if thresh(2) >= thresh(1)
  error('THRESH(2) must be smaller than THRESH(1)');
end

if ~isnumeric(d) || ~all(isfinite(d))
  error('D must be a numeric array of finite values');
end

if ~isstruct(dims) || (length(dims) ~= length(size(d)))
  error('DIMS must be a struct array that matches the dimensionality of D');
end
for i = 1:length(dims)
  if ~isnumeric(dims(i).grid) || ~isvector(dims(i).grid)
    error('DIMS must contain numeric grid vectors');
  end
  if any(diff(dims(i).grid) <= 0) || any(diff(diff(dims(i).grid)) ~= 0)
    error('grid spacings in DIMS must be uniform monotonic increasing');
  end
  if ~any(strcmp(dims(i).type,{'linear','circular'}))
    error('DIMS type must be either "linear" or "circular"');
  end
end

n = length(dims);
sz = size(d);
nvox = prod(sz);
[voxsubs{1:n}] = ind2sub(sz,(1:nvox)');
% build up a list of neighbor indices; each pixel faces 2*n neighbors (except
% for a few along the edges)
neighbors = zeros([nvox, 2*n]);
for i = 1:n
  tempsubs = voxsubs;
  if strcmp(dims(i).type,'linear')
    % edge voxels are assigned themselves as their own neighbors; these are
    % dummy entries which will be effective ignored in the flood-fill algorithm
    tempsubs{i} = voxsubs{i} - 1*(voxsubs{i} > 1);
    neighbors(:,2*i-1) = sub2ind(sz,tempsubs{:});
    tempsubs{i} = voxsubs{i} + 1*(voxsubs{i} < sz(i));
    neighbors(:,2*i) = sub2ind(sz,tempsubs{:});
  elseif strcmp(dims(i).type,'circular')
    % enforce wrap-around at edges
    tempsubs{i} = voxsubs{i} - 1;
    tempsubs{i}(tempsubs{i} == 0) = sz(i);
    neighbors(:,2*i-1) = sub2ind(sz,tempsubs{:});
    tempsubs{i} = voxsubs{i} + 1;
    tempsubs{i}(tempsubs{i} == sz(i)+1) = 1;
    neighbors(:,2*i) = sub2ind(sz,tempsubs{:});
  else
    error('dims type must be either "linear" or "circular"');
  end
end

% initialize output
h = zeros(sz);
% find max and min values in D
peak = max(d(:));
base = min(d(:));
% identify voxels that exceed thresh(1); we will use these 
% as starting seeds for the flood-fill algorithm
seeds = find(d(:) > thresh(1));
% running count of hot spots
count = 0;
while ~isempty(seeds)
  % initialize queue with the seed voxel that has the 
  % smallest subscript along the 1st dimension
  seedsubs{1} = ind2sub(sz,seeds);
  [junk, idx] = min(seedsubs{1});
  queue = seeds(idx);
  count = count+1;
  while ~isempty(queue)
    % mark all queued voxels as belonging to the current hotspot
    h(queue) = count;
    % replace queue with new members, which are neighbors of previously queued
    % voxels
    queue = neighbors(queue,:); 
    queue = unique(queue(:));
    % remove those voxels which are already marked
    queue(h(queue) > 0) = []; 
    % remove those voxels which do not exceed threshold
    queue(d(queue) <= thresh(2)) = [];
    % trim the list of seeds
    seeds(h(seeds) > 0) = [];
  end
end

