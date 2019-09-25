function [] = spikewaveforms(spikes, waves, timestamps,varargin)

%% plots waveforms for each channel, for a given subset of spikes (i.e.
%% spikes from a single epoch within a day)
% waves : all thresholded events, in xx-000 folder
% timestamps : timestamps of thresholded events, also in xx-000 folder
% spikes : timestamps of specific spikes you want


color=[.2 .2 .2];
if ~isempty(varargin)
    color=varargin{1};
end

spikes = spikes*10000;
waveforms = nan(40,4,length(spikes));           % initialize the collected set of waveforms

for s=1:length(spikes)                          % choose specified spikes (presumably from epoch) from all spikes from that tetrode & day
    i = find(timestamps==uint32(spikes(s)));             % retrieves index stored in timestamps
            waveforms(:,:,s) = waves(:,:,i);            % copies waveforms
end

seplot(1:40,(squeeze(waveforms(:,1,:)))','varType','std','color',color);

end



