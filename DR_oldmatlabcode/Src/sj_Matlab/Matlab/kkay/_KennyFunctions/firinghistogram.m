function [] = firinghistogram(spikes, pulses, windowsize, binsize)

% plots histogram given spike timestamps and stimulation pulse timestamps
% plots in non-overlapping bins of binsize ms

windowsize = windowsize/1000;     
binsize = binsize/1000;
nobins = windowsize/binsize+1;
binoffset = (nobins-1)/2;   %% to accommodate matlab's vector indexing
binnedout = zeros(nobins,1);   %% initializes output vector

avgfiringrate = length(spikes)/(pulses(length(pulses)) - pulses(1))

%% all arguments are now in timestamps of seconds


%% PSTH

for p=1:length(pulses)
    for b = -(nobins-1)/2:1:(nobins-1)/2                                                     %% iterates over all bins
        bincount = sum(spikes > (pulses(p)+b*binsize-binsize/2) &     ...                           %% bin counting of spikes
                       spikes <= (pulses(p)+b*binsize+binsize/2));
        binnedout(b+binoffset+1) = binnedout(b+binoffset+1) + bincount;                          %% adds running count
    end
end


t = -binsize*(nobins-1)/2:binsize:binsize*(nobins-1)/2;                                      %% time axis

baselinenobins = floor(0.25*length(binnedout));
baselinerate = sum(binnedout(1:baselinenobins))/ ...                                 %% baseline rate from first 0.25 of window
    (binsize*baselinenobins*length(pulses))                        


bar(t*1000,binnedout/((length(pulses)*binsize)*baselinerate));


%% spike waveform


spikewaveforms = nan(40,4,length(spikes));           % initialize

for s=1:length(spikes)                          % choose specified spikes (presumably from epoch) from all spikes from that tetrode & day
    find(waves(
    spikewaveforms = waves(:,:,s);      % copies waveforms
end

seplot(1:40,spikes(:,1,'varType','std'));

end

















