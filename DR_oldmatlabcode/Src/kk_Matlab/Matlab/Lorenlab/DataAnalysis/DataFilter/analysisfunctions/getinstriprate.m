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

  thresh = r.maxthresh(goodRips);

  p_starts = lookup(r.starttime(goodRips),p.data(:,1));
  p_mids = lookup(r.midtime(goodRips),p.data(:,1));
  p_ends = lookup(r.endtime(goodRips),p.data(:,1));

  rips = zeros(size(p.data,1),3);
  for i = 1:length(p_starts)
    rips(p_starts(i):p_ends(i),1) = thresh(i);
    rips(p_mids(i),2) = thresh(i);
  end

  noPosExcludePeriods = find(excludeperiods(:,1) > p.data(end,1) ...
    | excludeperiods(:,2) < p.data(1,1));
  p_starts = lookup(excludeperiods(~noPosExcludePeriods,1),p.data(:,1));
  p_ends = lookup(excludeperiods(~noPosExcludePeriods,2),p.data(:,1));
  for i = 1:length(p_starts)
    rips(p_starts(i):p_ends(i),:) = NaN;
  end

  % smooth rips
  if size(index,1) > 1
    out(t).rips = rips;
  else
    out = rips;
  end
end

