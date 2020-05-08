% xcorrstruct = spikexcorr(spiketimes1, spiketimes2, bin, tmax)
%    Computes the unnormalized cross correlation histogram of the spike
%    trains given by spiketimes1 and 2
%   Note that bin and tmax are given in seconds. tmax each side of center
%Returns a cross correlation structure with the following fields:
% time 	 	- the time correponding to the center of each bin  
% c1vsc2 	- the unnormalized cross-correlation histogram.
% nspikes1 	- the number of spikes from cell 1 
% nspikes2 	- the number of spikes from cell 2 

function [xcs] = spikexcorr(s1, s2, bin, tmax)

% adjust tmax upwards so that it is evenly divisible by the binsize
tmax = ceil(tmax / bin) * bin;

xcs.bin = bin;
xcs.tmax = tmax;
% do the cross correlation
xcs.nspikes1 = length(s1);
xcs.nspikes2 = length(s2);
% run the cross correlation.
if (xcs.nspikes1 & xcs.nspikes2)
    xcstmp = spikexcorrc(s1, s2, bin, tmax, 0);
    xcs.c1vsc2 = xcstmp.c1vsc2;
    xcs.time = (-tmax+bin/2):bin:(tmax-bin/2);
else
    xcs.c1vsc2 = [];
    xcs.time = [];
end
