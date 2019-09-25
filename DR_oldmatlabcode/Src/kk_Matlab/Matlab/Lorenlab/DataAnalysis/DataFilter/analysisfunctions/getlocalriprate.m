function [out] = getlocalriprate(index, times, excludeperiods, ripples, varargin)
% function [out] = getlocalriprate(index, times, excludeperiods, eeg, varargin)
%
%   index [day epoc tetrode]
%
%   out is windowed eeg organized as trial x time x tetrode

% assign the options

ratewindow = [-0.5 0];
[otherArgs] = procOptions(varargin);

out = [];
if isempty(times) | isempty(index)
  return;
end


goodPulses = find(~isExcluded(times, excludeperiods));
out = nan(length(goodPulses),size(index,1));

for t = 1:size(index,1) % number of tetrodes
  r = ripples{index(t,1)}{index(t,2)}{index(t,3)};
  for i = 1:length(goodPulses)
    out(i,t) = sum((r.starttime - times(goodPulses(i)) > ratewindow(1)) & ...
      (r.endtime - times(goodPulses(i)) < ratewindow(2)));
  end
end

