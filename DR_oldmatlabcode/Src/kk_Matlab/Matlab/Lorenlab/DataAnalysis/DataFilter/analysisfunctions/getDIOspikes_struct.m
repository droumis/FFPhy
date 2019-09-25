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

smoothing = 'none';
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

   if ~isempty(diofilter)
      filtresult = evaluatefilter2({pulses}, diofilter, 'activefield','pulseind');
      pulsetimes = pulses.pulsetimes(find(filtresult{1}(:,2))) / TIMESTAMPRATE;
      pulselength = pulses.pulselength(find(filtresult{1}(:,2)));
   else
      pulsetimes = pulses.pulsetimes(:,1) / TIMESTAMPRATE;
      pulselength = pulses.pulselength;
   end

   goodPulses = find(~isExcluded(pulsetimes, excludeperiods));

   % so now we have a list of pulse times.

   out(c).info = index(c,:);
   if isempty(goodPulses)
      if isempty(s) % this means there were no spikes from this neuron during this epoch
         out(c).psth = nan(0,L,1);
         out(c).meanpsth = nan(L,1);
         out(c).timestamps = {};
         out(c).pulsetimestamp = [];
         out(c).pulselength = [];
         out(c).theta = {};
      else
         out(c).psth = zeros(0,L,1);
         out(c).meanpsth = zeros(L,1);
         out(c).timestamps = {};
         out(c).pulsetimestamp = [];
         out(c).pulselength = [];
         out(c).theta = {};
      end
   else  % (if ~isempty(pulses.pulsetimes))
      out(c).psth = zeros(0,L,1);
      out(c).meanpsth = zeros(L,1);
      out(c).timestamps = {};
      out(c).pulsetimestamp = [];
      out(c).pulselength = [];
      out(c).theta = {};
      if ~isempty(s)
         out(c).psth = zeros(length(goodPulses),L,1);
         out(c).meanpsth = zeros(1,L);
         out(c).pulsetimestamp = pulsetimes(goodPulses);
         out(c).pulselength = pulselength(goodPulses);
         for i = 1:length(goodPulses)
            tmp_psth = zeros(L,1);
            tSt = (pulsetimes(goodPulses(i)) + pulsewin(1))*TIMESTAMPRATE;
            tEnd = (pulsetimes(goodPulses(i)) + pulsewin(end))*TIMESTAMPRATE;
            spkI = s.data(:,1)*TIMESTAMPRATE - tSt + 1; %(to account for 1-indexing)
            spkI = spkI * RESAMPRATE / TIMESTAMPRATE;
            spkI = round(spkI);

            inWindowInds = find((spkI > 0) & (spkI <= L));
            spkI = spkI(inWindowInds);
            tmp_psth(spkI) = 1;

            if ~isempty(spike_filter)
               tmp_psth = filter(spike_filter, 1, [tmp_psth; zeros(floor(length(spike_filter)/2),1)]);
               tmp_psth = tmp_psth(ceil(length(spike_filter)/2):end);
            end

            if ~isempty(spkI)
               out(c).psth(i,:) = tmp_psth;
               out(c).meanpsth = mean(out(c).psth);
               out(c).timestamps{i} = s.data(inWindowInds,1);
               if size(s.data,2) > 2
                  out(c).theta{i} = s.data(inWindowInds,3);
               end
            else
               out(c).timestamps{i} = s.data(inWindowInds,1);
               if size(s.data,2) > 2
                  out(c).theta{i} = s.data(inWindowInds,3);
               end
            end
         end
      end
   end
end


