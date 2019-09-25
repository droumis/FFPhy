function out = getpopulationevents(indices, excludeperiods, spikes, linpos, pos, ripples, cellinfo, cellcountthresh,varargin)
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
cellfilter = [];

for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})       
        switch(varargin{option})
            case 'cellfilter'
                cellfilter = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

index = [indices(1,1) indices(1,2)];
posdata = pos{index(1)}{index(2)}.data;
statematrix = linpos{index(1)}{index(2)}.statematrix;
%Find valid ripples
if isempty(cellfilter)
    riptimes = getripples(index, ripples, cellinfo, 'cellfilter', 'isequal($area, ''CA1'')','excludeperiods', excludeperiods,'minstd',3);
else
    riptimes = getripples(index, ripples, cellinfo, 'cellfilter', cellfilter,'excludeperiods', excludeperiods,'minstd',3);
end

out.index = [];
out.eventtraj = [];
out.eventdist = [];
out.eventtime = [];
out.eventimmobiletime = [];
out.eventdata = [];
out.peak = 0;
out.speed = [];
out.gamma = [];

spikecounts = [];
celldata = [];

if ~isempty(riptimes)
    %go through each cell and calculate the binned spike counts
    for cellcount = 1:size(indices,1)

        index = indices(cellcount,:);
        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
            if isfield(spikes{index(1)}{index(2)}{index(3)}{index(4)},'globalgammaphase')
                gam = spikes{index(1)}{index(2)}{index(3)}{index(4)}.globalgammaphase;
            else
                gam = nan(size(spiketimes));
            end
        else
            spiketimes = [];
            gam = [];
        end
        %Find valid spikes
        spiketimes = spiketimes(find(~isExcluded(spiketimes, excludeperiods)));
        spikebins = periodAssign(spiketimes, riptimes(:,[1 2]));
        
        if ~isempty(spiketimes)
            validspikes = find(spikebins);
            spiketimes = spiketimes(validspikes);
            spikebins = spikebins(validspikes);
            gam = gam(validspikes);
        end
        
        if ~isempty(spiketimes)
            tmpcelldata = [spiketimes spikebins gam];
            tmpcelldata(:,4) = cellcount;
        else
            tmpcelldata = [0 0 0 cellcount];
        end
        celldata = [celldata; tmpcelldata];
        spikecount = zeros(1,size(riptimes,1));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end

        spikecounts = [spikecounts; spikecount];
        out.index = [out.index; index];
    end
    %Sort all spikes by time
    celldata = sortrows(celldata,1);

    cellcounts = sum((spikecounts > 0));
    %Find all events with enough cells
    eventindex = find(cellcounts >= cellcountthresh);

    timeindex = lookup(riptimes(:,2),posdata(:,1));
    
    %Compute immobility time
    tmpvel = posdata(:,5);   
    tmpvel = tmpvel<(2);
    
    resetpoints = find(diff(tmpvel) < 0)+1;
    immobiletime = cumsum(tmpvel);
    for i = 1:length(resetpoints)
        immobiletime(resetpoints(i):end) = immobiletime(resetpoints(i):end)-immobiletime(resetpoints(i)-1);
    end
    immobiletime = immobiletime/30;
    immobile = immobiletime(timeindex);
      
    %For each event define: event time, ammount of time spent unmoving,
    %trajectory, linear distance, speed, and the times and identity of 
    %cells that fired during that event
    if iscell(statematrix.traj)
        traj = statematrix.traj{6}(timeindex);
    else
        traj = statematrix.traj(timeindex);
    end
    
    dist = statematrix.lindist(timeindex);
    vel = posdata(timeindex,8);

    for event = 1:length(eventindex)
        out.eventtime(event,1:2) = riptimes(eventindex(event),[1 2]);
        out.eventimmobiletime(event,1) = immobile(eventindex(event));
        out.eventtraj(event) = traj(eventindex(event));
        out.eventdist(event) = dist(eventindex(event));
        tmpind = find(celldata(:,2) == eventindex(event));
        out.eventdata(event).spiketimes = celldata(tmpind,1);
        out.eventdata(event).cellindex = celldata(tmpind,4);
        out.speed(event) = vel(eventindex(event));
        out.gamma(event).gammaphase = celldata(tmpind,3);
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



