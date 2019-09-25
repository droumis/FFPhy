function [out] = getclosestriptimes(index, times, excludeperiods, ripples, varargin)
% function [out] = getclosestriptimes(index, times, excludeperiods, ripples, varargin)
%
%   index [day epoc tetrode]
%
%   out is windowed eeg organized as trial x time x tetrode

% assign the options

pin = 48;
boundary = 0.01;
fields = {'closetimes','prerippower'};
preripwin = [-1 0];
[otherArgs] = procOptions(varargin);

out = [];
if isempty(times) | isempty(index)
  return;
end

goodPulses = find(~isExcluded(times, excludeperiods));

out = nan(length(goodPulses),length(fields),size(index,1));
for t = 1:size(index,1) % number of tetrodes
  rip = ripples{index(t,1)}{index(t,2)}{index(t,3)};
  k = 0;
  for i = 1:length(fields)
    k = k + 1;
    switch fields{i}
      case 'closetimes'
        inds = lookup(times(goodPulses),rip.endtime - 0.01,-1);
        out(:,k,t) = [rip.starttime(inds) - times(goodPulses)];
      case 'prerippower'
        for j = 1:length(goodPulses)
          prerips(j) = sum(rip.endtime > times(goodPulses(j)) + preripwin(1) & rip.endtime < times(goodPulses(j)) + preripwin(2));
        end
        out(:,k,t) = prerips; 
      otherwise
        error('field not supported');
      end
    end
  end
end

