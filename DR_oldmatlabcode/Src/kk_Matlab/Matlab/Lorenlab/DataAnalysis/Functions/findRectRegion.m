function [region, regionind] = findRectRegion(track_regions,pos)
% function [region, regionind] = findRectRegion(track_regions,pos)
%
% track_regions structure has two elements, 'name' and 'rect'
% 'rect' defines regions of track space, and 'name' is the
% corresponding name for each region.
%
% returns the name of the region for each [x,y] point in 'pos'
% and the corresponding index in the vector unique([track_regions.name])
% 
% NOTE: Overlapping regions will report only as the final one in
% the track_regions structure

possible_regions = unique([track_regions.name]);

region = nan(size(pos,1),1);
regionind = zeros(size(pos,1),1);

for j = 1:length(track_regions)
  rect = track_regions(j).rect;
  inside = ( pos(:,1) >= rect(1) ) & ...
      ( pos(:,1) <= rect(1) + rect(3) ) & ...
      ( pos(:,2) >= rect(2) ) & ...
      ( pos(:,2) <= rect(2) + rect(4) );
  region(inside,1) = track_regions(j).name;
  ind = find(possible_regions==track_regions(j).name);
  regionind(inside,1) = ind;
end

