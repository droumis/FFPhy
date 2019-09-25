function [out] = getDIOeeg(index, excludeperiods, spikes, DIO, varargin)
% function [out] = getDIOeeg(index, excludeperiods, spikes, DIO, varargin)
%
%   index [day epoc tetrode cell]
%
%   out is [DIO spike windows] or {DIO spike windows}

% assign the options

pulsewin = [-0.25 0.5];
pin = 48;
TIMESTAMPRATE = 10000;

RESAMPRATE = 1500;

smoothing = 'window';
smoothing_stddev = 30; % ms

spike_filter = [];

diofilter = '';
[otherArgs] = procOptions(varargin);

if strcmpi(smoothing,'gauss')
   spike_filter = exp(-0.5 * (linspace(-3,3,6*smoothing_stddev/1000*RESAMPRATE+1)).^2);
   spike_filter = spike_filter ./ sum(spike_filter);
elseif strcmpi(smoothing,'window')
   spike_filter = ones(10/1000*RESAMPRATE,1);
   spike_filter = spike_filter ./sum(spike_filter);
end

% convert pulsewin to timestamps
W = pulsewin * RESAMPRATE;
W = round(W(1)):round(W(2));
[zero_min,zero_ind] = find(min(abs(W)));
L = length(W);


for c = 1:size(index,1) % number of cell/tetrode/epochs
   s = spikes{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)};
   pulses = DIO{index(c,1)}{index(c,2)}{pin};

   if ~isempty(pulses.pulsetimes)
      if ~isempty(diofilter)
         filtresult = evaluatefilter2({pulses}, diofilter, 'activefield','pulseind');
         % we'll do everything in terms of the onset of the light
         pulsetimes = pulses.pulsetimes(find(filtresult{1}(:,2)));
      else
         pulsetimes = pulses.pulsetimes(:,1);
      end

      goodPulses = find(~isExcluded(pulsetimes, excludeperiods));

      % so now we have a list of pulse times.

      for i = 1:length(goodPulses)
         tmp_psth = zeros(L,1);
         tSt = pulsetimes(goodPulses(i)) + pulsewin(1)*TIMESTAMPRATE;
         tEnd = pulsetimes(goodPulses(i)) + pulsewin(end)*TIMESTAMPRATE;
         spkI = s.data(:,1)*TIMESTAMPRATE - tSt + 1; %(to account for 1-indexing)
         spkI = spkI * RESAMPRATE / TIMESTAMPRATE;
         spkI = round(spkI);
         spkI = spkI((spkI > 0) & (spkI <= L));
         tmp_psth(spkI) = 1;

         if ~isempty(spike_filter)
            tmp_psth = filter(spike_filter, 1, [tmp_psth; zeros(floor(length(spike_filter)/2),1)]);
            tmp_psth = tmp_psth(ceil(length(spike_filter)/2):end);
         end

         out(i,:,c) = tmp_psth;
      end
   else
      out(:,:,c) = zeros(0,L,1);
   end
end


