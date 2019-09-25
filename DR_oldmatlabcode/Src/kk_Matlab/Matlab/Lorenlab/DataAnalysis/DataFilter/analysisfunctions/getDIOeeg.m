function [out] = getDIOeeg(index, excludeperiods, eeg, DIO, varargin)
% out = getripdurations(index, excludeperiods, ripples, options)
%
%   index [day epoc tetrode]
%
%   out is [durations]

% assign the options

pulsewin = [-0.25 0.5];
pin = 48;
TIMESTAMPRATE = 10000;
diofilter = '';
[otherArgs] = procOptions(varargin);

out = [];

for t = 1:size(index,1) % number of tetrodes
   e = eeg{index(t,1)}{index(t,2)}{index(t,3)};
   pulses = DIO{index(t,1)}{index(t,2)}{pin};

   if ~isempty(pulses.pulsetimes)
      if ~isempty(diofilter)
         filtresult = evaluatefilter2({pulses}, diofilter, 'activefield','pulseind');
         pulsetimes = pulses.pulsetimes(find(filtresult{1}(:,2))) / TIMESTAMPRATE;
      else
         pulsetimes = pulses.pulsetimes(:,1) / TIMESTAMPRATE;
      end

      goodPulses = find(~isExcluded(pulsetimes, excludeperiods));

      etimes = e.starttime + [0:length(e.data)-1]*(1/e.samprate);

      startinds = lookup(pulsetimes(goodPulses),etimes);

      W = pulsewin * e.samprate;  % convert pulsewin from seconds to eeg samples
      W = round(W(1)):round(W(2));
      out(:,:,t) = nan(length(goodPulses),length(W),1);
      for i = 1:length(goodPulses)
         win = startinds(i) + W;
         if (win(1) >= 1) && (win(end) <= length(e.data))
            out(i,:,t) = e.data(win);
         elseif (win(end) > length(e.data))
            k = find(win > length(e.data),1);
            win = win(1:k-1);
            out(i,1:k-1,t) = e.data(win);
            % out(i,k+length(W),t) = nan;
         else
            k = find(win >= 1,1);
            win = win(k:end);
            % out(i,1:k-1,t) = nan;
            out(i,k:end,t) = e.data(win);
         end
      end
   % else
      % out(:,:,t) = [];
   end
end
if ~isempty(out)
   out = cat(2,out,repmat(pulsetimes(goodPulses),[1 1 size(index,1)]));
end

