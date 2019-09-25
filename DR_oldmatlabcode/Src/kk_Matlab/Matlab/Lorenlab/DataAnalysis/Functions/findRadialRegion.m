function [region, regionind] = findRadialRegion(track_regions,pos)
% function [region, regionind] = findRadialRegion(track_regions,pos)
%
% track_regions structure has elements, 'name' and radial data
% for each region
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
  p = [pos(:,1) - track_regions(j).center(1), ...
     pos(:,2) - track_regions(j).center(2)];
  a = atan2(p(:,2),p(:,1));
  r = sqrt(p(:,1).^2 +  p(:,2).^2);
  if diff(track_regions(j).arm_angles) < 0 % (angle rolling over)
     inside = (a >= track_regions(j).arm_angles(1)) | (a < track_regions(j).arm_angles(2));
  else
     inside = (a >= track_regions(j).arm_angles(1)) & (a < track_regions(j).arm_angles(2));
  end
  inside = inside & (r >= track_regions(j).rmin) & (r < track_regions(j).rmax);
  region(inside,1) = track_regions(j).name;

  ind = find(possible_regions == track_regions(j).name);
  regionind(inside,1) = ind;
end
