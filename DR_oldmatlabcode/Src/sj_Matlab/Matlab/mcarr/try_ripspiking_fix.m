%TRY TO FIX THIS.
load '/data21/mcarr/Bon/EEGnonreference/bonripple03-2-14.mat'
load '/data21/mcarr/Bon/bonspikes03.mat'
load '/data21/mcarr/Bon/bonripples03.mat'
load '/data21/mcarr/Bon/boncellinfo.mat'

day = 3; epoch = 2; tet = 1; c = 1;

eeg = ripple{day}{epoch}{tet};
spikes = spikes{day}{epoch}{tet}{c};

eegtimes = geteegtimes(eeg);
spikeind = lookup(spikes.data(:,1),eegtimes);
bin = 500;

%exclude spikes too near the edge
while spikeind(1) < bin
    spikeind(1) = [];
end
while spikeind(end) > length(eegtimes)-bin
    spikeind(end) = [];
end

%Pull out windows around each spike during ripples
riptimes = getripples([day epoch], ripples, cellinfo, 'cellfilter', ...
        '(isequal($area, ''CA1''))','minstd',3);
spikebins = periodAssign(spikes.data(:,1), riptimes(:,[1 2]));

validspikes = spikebins > 0;
xcorr_rip = zeros(sum(validspikes),2*bin+1);
spikeind = spikeind(validspikes);
for i = 1:length(spikeind)
    xcorr_rip(i,:) = ripple.data(spikeind(i)-bin:spikeind(i)+bin);
end

time = (median(diff(eegtimes))).*[-bin:1:bin];