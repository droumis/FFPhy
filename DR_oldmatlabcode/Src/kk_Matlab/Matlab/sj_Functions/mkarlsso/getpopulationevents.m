function out = getpopulationevents(indices, excludeperiods, spikes, linpos, window, cellcountthresh)
%
%Calculates the binned spike counts of all neurons in the index list.
%
%spikes - the 'spikes' cell array for the day you are analyzing
%linpos - the output of LINEARDAYPROCESS for the day you are analyzing. 
%indices - [day epoch tetrode cell]
%timebin- the length of each temporal bin (default 0.01 sec)
%excludeperiods - [start end] times for each exlcude period
%
%In the output, the spikecounts field is n by t, where n is the number of
%cells and t is the number of timebins.

if (nargin < 5)
    window = .3; %default window
end

index = [indices(1,1) indices(1,2)];
statematrix = linpos{index(1)}{index(2)}.statematrix;


timebins = min(statematrix.time):window:max(statematrix.time);

out.index = [];
out.eventtraj = [];
out.eventdist = [];
out.eventtime = [];
out.preeventcount = [];
out.eventdata = [];

spikecounts = [];
celldata = [];
%go through each cell and calculate the binned spike counts
for cellcount = 1:size(indices,1)

    index = indices(cellcount,:);
    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);      
    else
        spiketimes = [];
    end
    spiketimes = spiketimes(find(~isExcluded(spiketimes, excludeperiods)));
    spikebins = lookup(spiketimes,timebins);
    tmpcelldata = [spiketimes spikebins];
    tmpcelldata(:,3) = cellcount;
    celldata = [celldata; tmpcelldata];
    spikecount = zeros(1,length(timebins));
    for i = 1:length(spikebins)
        spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
    end
        
    spikecounts = [spikecounts; spikecount];
    out.index = [out.index; index];
end
celldata = sortrows(celldata,1); %sort all spikes by time




%newtimebinsInd = find(~isExcluded(timebins, excludeperiods));
%newtimes = timebins(newtimebinsInd);

cellcounts = sum((spikecounts > 0));
%cellcounts = cellcounts(newtimebinsInd);
eventindex = find(cellcounts >= cellcountthresh);

timeindex = lookup(timebins,statematrix.time);
traj = statematrix.traj(timeindex);
dist = statematrix.lindist(timeindex);


for event = 1:length(eventindex)
    out.eventtime(event) = timebins(eventindex(event));
    out.eventtraj(event) = traj(eventindex(event));
    out.eventdist(event) = dist(eventindex(event));
    try
        out.preeventcount(event) = cellcounts(eventindex(event)-1);
    catch
        out.preeventcount(event) = nan;
    end
    tmpind = find(celldata(:,2) == eventindex(event));
    out.eventdata(event).spiketimes = celldata(tmpind,1);
    out.eventdata(event).cellindex = celldata(tmpind,3);
end
    




