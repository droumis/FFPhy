function out = getradialrunregion(animaldir,animalprefix,epochs,varargin)
% function out = getradialrunregion(animaldir,animalprefix,epochs,varargin)
% out = getalllinstate(animaldir,animalprefix,epochs)
% Produces a cell structure with the fields:
% time, velocity (linearized)
%   EPOCHS - N by 2 matrix, columns are [day epoch]
%
%   OPTION: 'smooth', default no smoothing
%                   to compute linear speed, can smooth linear position data.
%                   It is smoothed with a gaussian of length VSW and std VSW/4.
%                   default for lineardayprocess is 2 seconds

% smooth = [];
% if ~isempty(varargin)
%     smooth = varargin{2};
% end

loaddays = unique(epochs(:,1));
pos = loaddatastruct(animaldir, animalprefix, 'pos', loaddays);
load(fullfile(animaldir,[animalprefix,'track_regions']));

for i = 1:size(epochs,1)
   out{epochs(i,1)}{epochs(i,2)}.time = pos{epochs(i,1)}{epochs(i,2)}.data(:,1);
   out{epochs(i,1)}{epochs(i,2)}.region = repmat(' ',length(out{epochs(i,1)}{epochs(i,2)}.time),1);
   for j = 1:length(track_regions)
      p = [pos{epochs(i,1)}{epochs(i,2)}.data(:,2) - track_regions(j).center(1), ...
         pos{epochs(i,1)}{epochs(i,2)}.data(:,3) - track_regions(j).center(2)];
      a = atan2(p(:,2),p(:,1));
      r = sqrt(p(:,1).^2 +  p(:,2).^2);
      if diff(track_regions(j).arm_angles) < 0 % (angle rolling over)
         inside = (a >= track_regions(j).arm_angles(1)) | (a < track_regions(j).arm_angles(2));
      else
         inside = (a >= track_regions(j).arm_angles(1)) & (a < track_regions(j).arm_angles(2));
      end
      inside = inside & (r >= track_regions(j).rmin) & (r < track_regions(j).rmax);
      out{epochs(i,1)}{epochs(i,2)}.region(inside,1) = track_regions(j).name;
    end
end

end
