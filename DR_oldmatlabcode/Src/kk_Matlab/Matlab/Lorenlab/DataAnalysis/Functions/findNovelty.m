function out = findNovelty(track_regions, task, position)
% function out = findNovelty(animalinfo, task, position)
%
% animalinfo.regionfunction should specify what function to call
% to determine what region of the track is novel.
%
% task.exposure is a 2x1 vector: 
%   [day of novel arm exposure, run of novel arm exposure]
%   [2 1] = day 2 on novel arm, run 1
%   [1 2] = day 1 on novel arm, run 2
%   [inf 3] = all arms familiar, run 3 in this configuration today

if isfield(track_regions,'rect') % rectangular track
  [region,regionind] = findRectRegion(track_regions,position);
elseif isfield(track_regions,'arm_angles') % radial track
  [region,regionind] = findRadialRegion(track_regions,position);
else
  error('findNovelty: Unknown track region structure.');
end

if ~isfield(task,'geo_novelty') 
  out = nan(size(regionind,1),2); 
  % warning('Task structure does not have expected geo_novelty field');
  return;
end

out = inf(size(regionind,1),2); % default is FAMILIAR
for i = 1:length(regionind)
  out(i,:) = task.geo_novelty(regionind(i),:);
end

