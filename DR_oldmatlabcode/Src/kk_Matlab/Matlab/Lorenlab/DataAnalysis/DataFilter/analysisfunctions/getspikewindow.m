function [out] = getspikewindow(index, times, excludeperiods, spikes, varargin)
% function [out] = getspikewindow(index, times, excludeperiods, spikes, varargin)
%
%   index [day epoc tetrode cell]
%
%   out is windowed eeg organized as trial x time x tetrode

% spikewin = [-0.25 0.5]; % size of window around DIO pulses
% TIMESTAMPRATE = 10000;
% RESAMPRATE = 1500; - e.g., EEG sampling rate
% smoothing = 'window'; % also 'gauss', '' = no smoothing
% smoothing_stddev = 2; % ms
% output = 'psth'; % could be 'spikes'
% inputType = 'spike'; % could be 'multi'
% spike_filter = []; % used instead of default smoothing

% assign the options

spikewin = [-0.25 0.5]; % size of window around DIO pulses

TIMESTAMPRATE = 10000;
RESAMPRATE = 1500;

smoothing = 'window';
smoothing_stddev = 2; % ms

output = 'psth';

inputType = 'spike';

spike_filter = [];
[otherArgs] = procOptions(varargin);

out = [];
if isempty(times) | isempty(index)
  return;
end

if ~(strcmpi(output,'psth') | strcmpi(output,'spikes'))
  error('Unknown output type for getspikewindow.');
end

if isempty(spike_filter)
  if strcmpi(smoothing,'gauss')
     spike_filter = exp(-0.5 * (linspace(-3,3,6*smoothing_stddev/1000*RESAMPRATE+1)).^2);
     spike_filter = spike_filter ./ sum(spike_filter);
  elseif strcmpi(smoothing,'window')
     spike_filter = ones(10/1000*RESAMPRATE,1);
     spike_filter = spike_filter ./sum(spike_filter);
  end
end

% convert spikewin to timestamps
% W = spikewin * RESAMPRATE;
% W = round(W(1)):round(W(2));
% [zero_min,zero_ind] = find(min(abs(W)));
% L = length(W);
L = ceil((spikewin(2)-spikewin(1)) * RESAMPRATE);

for c = 1:size(index,1) % number of cell/tetrode/epochs
  if strcmp(inputType,'spike')
    s = spikes{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}.data;
  elseif strcmp(inputType, 'multi')
    s = spikes{index(c,1)}{index(c,2)}{index(c,3)}/TIMESTAMPRATE;
  end

  if isempty(s)
    continue;
  end

  goodPulses = find(~isExcluded(times, excludeperiods));

  if strcmpi(output,'psth')
    out(:,:,c) = nan(length(goodPulses),L,1);
  end
  for i = 1:length(goodPulses)

    % REPLACE WITH SPIKEXCORR!!!!!


     tSt = times(goodPulses(i)) + spikewin(1);
     tEnd = times(goodPulses(i)) + spikewin(end);
     if strcmpi(output,'spikes')
       out(c).timestamps{i} = s(s(:,1)>tSt & s(:,1)<=tEnd,1) - times(goodPulses(i));
       out(c).pulsetime(i) = times(goodPulses(i));
       continue;
     end

     % NOTE - THIS APPROACH IS CONSERVATIVE IN THE FOLLOWING WAY:
     %  SPIKES IN BIN 1 ARE GUARANTEED TO HAVE OCCURED FOLLOWING
     %  tSt, SPIKES IN BIN 2, 1/RESAMPRATE s following tSt,
     %  etc.
     tmp_psth = zeros(L,1);
     spkI = s(:,1) - tSt;
     spkI = spkI * RESAMPRATE + 1; %(to account for 1-indexing)
     spkI = round(spkI);
     spkI = spkI((spkI > 0) & (spkI <= L));
     tmp_psth(spkI) = 1;

     % NON CAUSAL
     if ~isempty(spike_filter)
        tmp_psth = filter(spike_filter, 1, [tmp_psth; zeros(floor(length(spike_filter)/2),1)]);
        tmp_psth = tmp_psth(ceil(length(spike_filter)/2):end);
     end

     if strcmpi(output,'psth')
       out(i,:,c) = tmp_psth;
     end
  end

end

