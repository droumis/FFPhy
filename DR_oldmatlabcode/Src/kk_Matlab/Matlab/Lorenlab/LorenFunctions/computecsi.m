function [csi propb] = computecsi(s, amp, blen) 
% function [csi propburst] = computecsi(spiketimes, spike_amplitudes, burstlen) 
%	Computes 
%		csi  - the complex spike index and
%		propburst - the proportion of spikes that are part of bursts
%
% 	based on the inputed set of times, amplitudes and the burst length 
%	in milliseconds.

% check that the burst length is reasonable
if (blen < 1)
    error('burst length must be > 1 ms');
end

csi = 0;
propb = 0;

% convert blen to seconds
blen = blen / 1000;

nspikes = length(s);
if (nspikes < 2)
    return;
end

%get the  list of intervals;
tmpisi = diff(s);

% find the intervals less than the burst length
bspikes = find(tmpisi < blen);

% for each pair burst spikes, determine whether the amplitude of the spikes
% decreased
ampdiff = amp(bspikes) - amp(bspikes+1);

% the csi is the number of ampdiff that are greater than zero divided by the
% total number of spikes
csi = length(find(ampdiff > 0)) / nspikes;

% a single interval corresponds to two spikes, but if there are two or more
% intervals next to each other, subtract one for each adjacent interaval
nadjacent = 0;
if (length(bspikes) > 2)
        nadjacent = length(find(diff(bspikes) == 1));
end
propb = (length(bspikes) * 2 - nadjacent) / nspikes;
