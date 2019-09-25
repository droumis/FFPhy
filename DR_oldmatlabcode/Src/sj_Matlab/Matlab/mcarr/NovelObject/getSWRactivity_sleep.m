function out = getSWRactivity_sleep(indices, excludeperiods, spikes, pos, ripples, cellinfo ,varargin)

%Calculates activation probability and the mean rate during SWRs

%spikes - the 'spikes' cell array for the day you are analyzing
%indices - [day epoch tetrode cell]
%timebin- the length of each temporal bin (default 0.01 sec)
%excludeperiods - [start end] times for each exlcude period
%

out.index = [];
out.totalSWRs = [];
out.totalSWRduration = [];
out.activationprobability = [];
out.SWRrate = [];
out.totalduration = [];

cellfilter = [];
for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'cellfilter'
                cellfilter = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    end
end

index = [indices(1,1) indices(1,2)];
pos = pos{index(1)}{index(2)};

%Find valid ripples
if ~isempty(evaluatefilter(cellinfo{index(1)}{index(2)},'isequal($area,''CA1'')'))
    if isempty(cellfilter)
        riptimes = getripples(index, ripples, cellinfo, 'cellfilter', 'isequal($area, ''CA1'')','excludeperiods', excludeperiods,'minstd',3);
    else
        riptimes = getripples(index, ripples, cellinfo, 'cellfilter', cellfilter, 'excludeperiods',excludeperiods, 'minstd', 3);
    end

    out.totalSWRs = size(riptimes,1); out.totalSWRduration = sum(riptimes(:,2)-riptimes(:,1));
    out.totalduration = pos.data(end,1)-pos.data(1,1) - sum(excludeperiods(:,2)-excludeperiods(:,1));

    %Try limiting to first five minutes for all sessions
    fivemin = pos.data(1,1) + 5*60;
    riptimes(riptimes(:,3)>fivemin,:) = [];
    out.totalSWRs = size(riptimes,1); out.totalSWRduration = sum(riptimes(:,2)-riptimes(:,1));
    out.totalduration = 5*60 - sum(excludeperiods(:,2)-excludeperiods(:,1));

    if ~isempty(riptimes)
        %go through each cell
        for cellcount = 1:size(indices,1)

            index = indices(cellcount,:);
            if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
                spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
            else
                spiketimes = [];
            end

            spiketimes = spiketimes(find(~isExcluded(spiketimes, excludeperiods)));
            spikebins = periodAssign(spiketimes, riptimes(:,[1 2]));

            if ~isempty(spiketimes)
                validspikes = find(spikebins);
                spikebins = spikebins(validspikes);
            end

            spikecount = zeros(1,size(riptimes,1));
            for i = 1:length(spikebins)
                spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
            end

            %Look at activation probability and SWR rate
            out.index = [out.index; index];
            out.activationprobability(cellcount) = sum(spikecount>0)/out.totalSWRs;
            out.SWRrate(cellcount) = sum(spikecount)/out.totalSWRs;
        end    
    end
end