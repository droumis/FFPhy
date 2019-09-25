function out = getpopulationrates(indices, excludeperiods, spikes, linpos, timebin)
%trajdata = FILTERCALCLINFIELDS(index, excludeperiods, spikes, linpos)
%trajdata = FILTERCALCLINFIELDS(index, excludeperiods, spikes, linpos, binsize)
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
    timebin = .01; %default time bin
end

index = [indices(1,1) indices(1,2)];
statematrix = linpos{index(1)}{index(2)}.statematrix;


tmptimebins = min(statematrix.time):timebin:max(statematrix.time);
out.timebins = tmptimebins(find(~isExcluded(tmptimebins, excludeperiods)));
timeindex = lookup(out.timebins,statematrix.time);
out.traj = statematrix.traj(timeindex);
out.dist = statematrix.lindist(timeindex);
out.index = [];
out.spikecounts = uint8([]);

%go through each cell and calculate the binned spike counts
for cellcount = 1:size(indices,1)

    index = indices(cellcount,:);
    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
        spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);      
    else
        spiketimes = [];
    end
    spiketimes = spiketimes(find(~isExcluded(spiketimes, excludeperiods)));
    spikebins = lookup(spiketimes,out.timebins);
    spikecount = zeros(1,length(out.timebins));
    for i = 1:length(spikebins)
        spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
    end
        
    out.spikecounts = [out.spikecounts;uint8(spikecount)];
    out.index = [out.index; index];
end







