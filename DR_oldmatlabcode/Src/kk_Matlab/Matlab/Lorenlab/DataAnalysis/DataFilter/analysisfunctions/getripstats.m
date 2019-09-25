function [out] = getripplestats(index, times, excludeperiods, ripple, ripples, varargin)
% function [out] = getripplestats(index, excludeperiods, ripple, varargin)
%
%   index [day epoc tetrode]
%
%   out is windowed ripple organized as trial x time x tetrode

% assign the options

smoothing_width = 0.004;

[otherArgs] = procOptions(varargin);

out = [];
if isempty(times) | isempty(index)
  return;
end

goodPulses = find(~isExcluded(times(:,1), excludeperiods));

for t = 1:size(index,1) % number of tetrodes
  r = ripple{index(t,1)}{index(t,2)}{index(t,3)};

  kernel = gaussian(smoothing_width*r.samprate, ceil(8*smoothing_width*r.samprate));
  renv = smoothvect(double(r.data(:,2)), kernel);
  env_mean = mean(renv);
  env_std = std(abs(renv));

  rtimes = r.starttime + [0:length(renv)-1]*(1/r.samprate);

  startinds = lookup(times(goodPulses,1),rtimes);
  endinds = lookup(times(goodPulses,2),rtimes);

  out(:,:,t) = nan(length(goodPulses),3,1);
  for i = 1:length(goodPulses)
    win = startinds(i):endinds(i);
    % win = (startinds(i)-round(50*1.5)):(endinds(i)+round(50*1.5));

    out(i,1,t) = (max(renv(win)))/env_std; % peak
    out(i,2,t) = (mean(renv(win)))/env_std; % average
    out(i,3,t) = (sum(renv(win)))/env_std; % total
  end

end

