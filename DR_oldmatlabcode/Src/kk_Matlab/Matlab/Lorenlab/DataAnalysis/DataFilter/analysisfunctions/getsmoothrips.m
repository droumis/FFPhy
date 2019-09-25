function [out] = getsmoothrips(index, excludeperiods, ripples, pos, varargin)
% function [out] = getsmoothrips(index, excludeperiods, ripples, pos, varargin)
%
%   index [day epoc tetrode]
%
%   out is struct array of ripples

% assign the options

minthresh = 2;
smoothRate = 60; % corresponds to roughly 2 seconds

[otherArgs] = procOptions(varargin);

out = [];
if isempty(index)
  return;
end

if size(unique(index(:,1:2),'rows'),1) ~= 1
  error('getsmoothriprate only takes one epoch at a time');
end

p = pos{index(1,1)}{index(1,2)};


ds_t = downsample(p.data(:,1),smoothRate,floor(smoothRate/2));
ds_v = filter(ones(smoothRate,1)/smoothRate,1,[zeros(floor(smoothRate/2),1);p.data(:,5)]);
ds_v = ds_v(floor(smoothRate/2)+1:end);
ds_v = downsample(ds_v,smoothRate,floor(smoothRate/2));
% ds_t = decimate(p.data(:,1),smoothRate);
% ds_v = decimate(p.data(:,5),smoothRate);

inds = find(ds_v > 100);
if ~isempty(inds)
  keyboard
end


for t = 1:size(index,1) % number of tetrodes
  r = ripples{index(t,1)}{index(t,2)}{index(t,3)};

  aboveThresh = r.maxthresh > minthresh;

  excludedStarts = isExcluded(r.starttime, excludeperiods);
  excludedEnds = isExcluded(r.endtime, excludeperiods);
  excludedNoPos = r.starttime < p.data(1,1) | r.endtime > p.data(end,1);
  goodRips = find(~(excludedStarts | excludedEnds | excludedNoPos) & aboveThresh);

  out(:,1,t) = hist(r.midtime(goodRips),ds_t);
  % out(:,1,t) = zeros(length(ds_t),1);
  % t_inds = lookup(r.midtime(goodRips),ds_t);
  % for i = 1:length(t_inds)
    % out(t_ind(i),1,t) = out(t_ind(i),1,t) + 1;
  % end

  out(:,2,t) = ds_v;

end

