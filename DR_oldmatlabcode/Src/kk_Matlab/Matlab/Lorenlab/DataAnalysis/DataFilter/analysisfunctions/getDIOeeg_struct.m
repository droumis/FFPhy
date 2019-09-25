function [out] = getDIOeeg(index, excludeperiods, eeg, DIO, varargin)
% out = getripdurations(index, excludeperiods, ripples, options)
%
%   index [day epoc tetrode]
%
%   out is [durations]

% assign the options

if length(varargin) > 0
   if iscell(varargin{1})
      ripples = varargin{1};
      varargin = {varargin{2:end}};
   else
      ripples = [];
   end
end

pulsewin = [-0.25 0.5];
pin = 48;
TIMESTAMPRATE = 10000;
diofilter = '';
closestripples = 0;
[otherArgs] = procOptions(varargin);


if closestripples ~= 0 & isempty(ripples)
   error('To use closestripples option, must pass ripples struct to function.');
end

for t = 1:size(index,1) % number of day/epoch/tetrodes
%    e = eeg{index(t,1)}{index(t,2)}{index(t,3)};
   S = eeg{index(t,1)}{index(t,2)}{index(t,3)}.samprate;
   L = length(eeg{index(t,1)}{index(t,2)}{index(t,3)}.data);
   pulses = DIO{index(t,1)}{index(t,2)}{pin};

   W = pulsewin * S;  % convert pulsewin from seconds to eeg samples
   W = round(W(1)):round(W(2));


   if ~isempty(pulses.pulsetimes)
      if ~isempty(diofilter)
         filtresult = evaluatefilter2({pulses}, diofilter, 'activefield','pulseind');
         filteredInds = find(filtresult{1}(:,2));
         pulsetimes = pulses.pulsetimes(filteredInds) / TIMESTAMPRATE;
      else
         pulsetimes = pulses.pulsetimes(:,1) / TIMESTAMPRATE;
         filteredInds = [1:length(pulsetimes)];
      end

      goodPulses = find(~isExcluded(pulsetimes, excludeperiods));

      if ~isempty(goodPulses)
         etimes = eeg{index(t,1)}{index(t,2)}{index(t,3)}.starttime + [0:L-1]*(1/S);
%          etimes = e.starttime + [0:length(e.data)-1]*(1/e.samprate);

         startinds = lookup(pulsetimes(goodPulses),etimes);
         clear etimes

         out.data(:,:,t) = nan(length(goodPulses),length(W),1);
         out.mindata = nan(length(goodPulses),1);
         out.pulsedata = [pulses.pulselength(filteredInds) pulses.frequency(filteredInds) pulses.frequency2(filteredInds) pulses.index_in_sequence(filteredInds) pulses.timesincelast(filteredInds)];
         if closestripples
            out.closestripple(:,:,t) = nan(length(goodPulses),3);
         end
         for i = 1:length(goodPulses)
            win = startinds(i) + W;
            out.inds(i,:) = [win(1) win(end)];
            out.times(i,:) = out.inds(i,:) / S;
            if (win(1) >= 1) && (win(end) <= L)
               out.data(i,:,t) = eeg{index(t,1)}{index(t,2)}{index(t,3)}.data(win);
            elseif (win(end) > L)
               k = find(win > L,1);
               win = win(1:k-1);
               out.data(i,1:k-1,t) = eeg{index(t,1)}{index(t,2)}{index(t,3)}.data(win);
            else
               k = find(win >= 1,1);
               win = win(k:end);
               out.data(i,k:end,t) = eeg{index(t,1)}{index(t,2)}{index(t,3)}.data(win);
            end
            out.mindata(i) = min(out.data(i,:));

            if closestripples
               tpulse = pulsetimes(goodPulses(i));
               closest_rip = find(ripples{index(t,1)}{index(t,2)}{index(t,3)}.starttime > tpulse,1);
               if isempty(closest_rip)
                  out.closestripple(i,:,t) = [nan nan nan];
               elseif closest_rip < closestripples
                  out.closestripple(i,:,t) = [nan nan nan];
               else
                  ii = closest_rip - closestripples;
                  out.closestripple(i,:,t) = [
                     ripples{index(t,1)}{index(t,2)}{index(t,3)}.starttime(ii)-tpulse  ...
                     ripples{index(t,1)}{index(t,2)}{index(t,3)}.peak(ii) ...
                     ripples{index(t,1)}{index(t,2)}{index(t,3)}.energy(ii)];
               end
            end
         end
      else
         out.data(:,:,t) = nan(0,length(W),1);
         out.mindata = nan(0,1);
         out.pulsedata = nan(0,5);
         if closestripples
            out.closestripples(:,:,t) = nan(0,3);
         end
         out.inds = nan(0,2);
         out.times = nan(0,2);
      end
   else
      out.data(:,:,t) = nan(0,length(W),1);
      out.mindata = nan(0,1);
      out.pulsedata = nan(0,5);
      if closestripples
         out.closestripples(:,:,t)= nan(0,3);
      end
      out.inds = nan(0,2);
      out.times = nan(0,2);
   end

   eeg{index(t,1)}{index(t,2)}{index(t,3)}.data = []; % clear up memory as we go
end
if ~isempty(out.data)
   out.pulsetimes = pulsetimes(goodPulses);
else
   out.pulsetimes = [];
end

