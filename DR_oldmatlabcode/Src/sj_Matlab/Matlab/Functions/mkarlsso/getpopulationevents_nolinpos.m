function out = getpopulationevents_nolinpos(indices, excludeperiods, spikes, pos, ripples, cellinfo, cellcountthresh)
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


index = [indices(1,1) indices(1,2)];
posdata = pos{index(1)}{index(2)}.data;
[riptimes] = getripples(index, ripples, cellinfo, 'cellfilter', '(isequal($area, ''CA1'')|isequal($area, ''CA3''))','excludeperiods', excludeperiods,'minstd',3);
%[riptimes] = getripples(index, ripples, cellinfo, 'cellfilter', '(isequal($area, ''CA1''))','excludeperiods', excludeperiods,'minstd',3);

%riptimes(:,1) = riptimes(:,2)-(window/2);
%riptimes(:,3) = riptimes(:,2)+(window/2);


out.index = [];
out.eventtraj = [];
out.eventdist = [];
out.eventtime = [];
out.eventimmobiletime = [];
out.std = [];
out.preeventcount = [];
out.eventdata = [];
out.peak = 0;

spikecounts = [];
celldata = [];

if ~isempty(riptimes)
    %go through each cell and calculate the binned spike counts
    for cellcount = 1:size(indices,1)

        index = indices(cellcount,:);
        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
        else
            spiketimes = [];
        end
        spiketimes = spiketimes(find(~isExcluded(spiketimes, excludeperiods)));
        spikebins = periodAssign(spiketimes, riptimes(:,[1 2]));
        
        %spikebins = lookup(spiketimes,riptimes(:,2));
        if ~isempty(spiketimes)
            validspikes = find(spikebins);
            spiketimes = spiketimes(validspikes);
            spikebins = spikebins(validspikes);
            
%             spikedeviation = abs(spiketimes - riptimes(spikebins,2));
%             validspikes = find(spikedeviation <= (window/2));
%             spiketimes = spiketimes(validspikes);
%             spikebins = spikebins(validspikes);
        end
        tmpcelldata = [spiketimes spikebins];
        tmpcelldata(:,3) = cellcount;
        celldata = [celldata; tmpcelldata];
        spikecount = zeros(1,size(riptimes,1));
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

    timeindex = lookup(riptimes(:,2),posdata(:,1));
   
    tmpvel = posdata(:,5);
    %tmpvel = tmpvel<(max(tmpvel)*.05);
    tmpvel = tmpvel<(2);
    
 
    
    resetpoints = find(diff(tmpvel) < 0)+1;
    immobiletime = cumsum(tmpvel);
    for i = 1:length(resetpoints)
        immobiletime(resetpoints(i):end) = immobiletime(resetpoints(i):end)-immobiletime(resetpoints(i)-1);
    end
    immobiletime = immobiletime/30;
    immobile = immobiletime(timeindex);
    


    for event = 1:length(eventindex)
        out.eventtime(event,1:2) = riptimes(eventindex(event),[1 2]);
        out.eventimmobiletime(event,1) = immobile(eventindex(event));
        tmpind = find(celldata(:,2) == eventindex(event));
        out.eventdata(event).spiketimes = celldata(tmpind,1);
        out.eventdata(event).cellindex = celldata(tmpind,3);
        %out.std(event,1) = ripstd(eventindex(event));
    end
else
    warning('No ripples found')
end
    


function out = periodAssign(times, periods)
% out = periodAssign(times, periods)
% TIMES is a vector of times
% PERIODS is an N by 2 list of start and end times
% Returns the index of the period that each time falls into.  If a time
% does not fall into one of the periods, a zero is returned.
% This function assumes that the periods are not overlapping.
%

if ~isempty(periods)
    oneborder = [(periods(:,1)-.0000001);periods(:,2)+.0000001];
    oneborder(:,2) = 0;
    insideborder = [(periods(:,1)+.0000001) (1:length(periods))'; (periods(:,2)-.0000001) (1:length(periods))'];    
    sortedMatrix = [[-inf 0]; sortrows([oneborder;insideborder],1); [inf 0]];
else
    sortedMatrix = [[-inf 0]; [inf 0]];
end
out = sortedMatrix(lookup(times,sortedMatrix(:,1)),2);
