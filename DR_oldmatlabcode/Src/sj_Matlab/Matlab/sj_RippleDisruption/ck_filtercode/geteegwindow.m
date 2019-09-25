function [out] = geteegwindow(index, times, excludeperiods, eeg, varargin)
% function [out] = geteegwindow(index, excludeperiods, eeg, varargin)
%
%   index [day epoc tetrode]
%
%   out is windowed eeg organized as trial x time x tetrode

% assign the options

eegwindow = [-0.25 0.5];
[otherArgs] = procOptions(varargin);

out = [];
if isempty(times) | isempty(index)
  return;
end


for t = 1:size(index,1) % number of tetrodes
  e = eeg{index(t,1)}{index(t,2)}{index(t,3)};

  if (size(e.data,2) > 1) % (complex - phase / magnitude)
    d = exp(j*e.data(:,1)) .* e.data(:,2);
  else
    d = e.data;
  end

  goodPulses = ~isExcluded(times(:,1), excludeperiods);
  if size(times,2) > 1
    for tt = 1:size(times,2)
      goodPulses = goodPulses & ~isExcluded(times(:,1), excludeperiods);
    end
  end
  goodPulses = find(goodPulses);


  etimes = e.starttime + [0:length(d)-1]*(1/e.samprate);

  startinds = lookup(times(goodPulses,1),etimes);

  W = eegwindow * round(e.samprate);  % convert eegwindow from seconds to eeg samples
  W = round(W(1)):round(W(2));

  out(:,:,t) = nan(length(goodPulses),length(W),1);
  for i = 1:length(goodPulses)
    win = startinds(i) + W;
    if (win(1) >= 1) && (win(end) <= length(d))
      out(i,:,t) = d(win);
    elseif (win(end) > length(d))
      k = find(win > length(d),1);
      win = win(1:k-1);
      out(i,1:k-1,t) = d(win);
    else
      k = find(win >= 1,1);
      win = win(k:end);
      % out(i,1:k-1,t) = nan;
      out(i,k:end,t) = d(win);
    end
  end
end

