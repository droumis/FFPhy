function out = getrunregion(animaldir,animalprefix,epochs,varargin)
% function out = getrunregion(animaldir,animalprefix,epochs,varargin)
% Produces a cell structure with the fields:
% time, velocity (linearized)
%   EPOCHS - N by 2 matrix, columns are [day epoch]
%

loaddays = unique(epochs(:,1));
pos = loaddatastruct(animaldir, animalprefix, 'pos', loaddays);
load(fullfile(animaldir,[animalprefix,'track_regions']));
track_rects = vertcat(track_regions(:).rect);

% rect is [left bottom width height]

for i = 1:size(epochs,1)
   out{epochs(i,1)}{epochs(i,2)}.time = pos{epochs(i,1)}{epochs(i,2)}.data(:,1);
   out{epochs(i,1)}{epochs(i,2)}.region = repmat(' ',length(out{epochs(i,1)}{epochs(i,2)}.time),1);
   for j = 1:size(track_rects,1)
      inside = ( pos{epochs(i,1)}{epochs(i,2)}.data(:,2) >= track_rects(j,1) ) & ...
        ( pos{epochs(i,1)}{epochs(i,2)}.data(:,2) <= track_rects(j,1) + track_rects(j,3) ) & ...
        ( pos{epochs(i,1)}{epochs(i,2)}.data(:,3) >= track_rects(j,2) ) & ...
        ( pos{epochs(i,1)}{epochs(i,2)}.data(:,3) <= track_rects(j,2) + track_rects(j,4) );
      out{epochs(i,1)}{epochs(i,2)}.region(inside,1) = track_regions(j).name;
    end
end
