function [out] = getrippos(index, excludeperiods, ripples, pos, varargin)
% function [out] = getrippos(index, excludeperiods, ripples, pos, varargin)
%
%   index [day epoc tetrode]
%
%   out is struct array of ripples

% assign the options

minthresh = 2;

[otherArgs] = procOptions(varargin);

out = [];
if isempty(index)
  return;
end

for t = 1:size(index,1) % number of tetrodes
  r = ripples{index(t,1)}{index(t,2)}{index(t,3)};

  p = pos{index(t,1)}{index(t,2)};

  aboveThresh = r.maxthresh > minthresh;

  excludedStarts = isExcluded(r.starttime, excludeperiods);
  excludedEnds = isExcluded(r.endtime, excludeperiods);
  excludedNoPos = r.starttime < p.data(1,1) | r.endtime > p.data(end,1);
  goodRips = find(~(excludedStarts | excludedEnds | excludedNoPos) & aboveThresh);

  p_inds = lookup(r.midtime(goodRips),p.data(:,1));
  out(t).pos = p.data(p_inds,:);

  out(t).starttime = r.starttime(goodRips);
  out(t).midtime = r.midtime(goodRips);
  out(t).peaktime = r.peaktime(goodRips);
  out(t).endtime = r.endtime(goodRips);
  out(t).duration = r.endtime(goodRips) - r.starttime(goodRips);
  out(t).peak = r.peak(goodRips);
  out(t).energy = r.energy(goodRips);
  out(t).maxthresh = r.maxthresh(goodRips);
end

