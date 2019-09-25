function out = calcripplespikingmodulation(index, excludeperiods, spikes, ripples,cellinfo, varargin)
% Computes the modulation of spikes by beta, gamma, and ripple oscillations during SWRs.
% out.celldata is a Nx6 matrix with N rows for each spike and 6 rows
%   1: spiketime
%   2: which ripple the spike occured in
%   3: max ripple thresh
%   4: gamma phase
%   5: ripple phase
%   6: which cell
%   7: area recorded

%Set options
cellfilter = [];
min_cells = 1;

celldata = [];
spikecounts = [];
minthresh = 3;

%Process options
for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})
        switch(varargin{option})
            case 'cellfilter'
                cellfilter = varargin{option+1};
            case 'min_cells'
                min_cells = varargin{option+1};
            case 'minthresh'
                minthresh = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

%Find valid ripples
if isempty(cellfilter)
    [riptimes ripstd]= getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
        '(isequal($area, ''CA1''))',...
        'excludeperiods', excludeperiods,'minstd',minthresh);
else
    [riptimes ripstd]= getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
        cellfilter,'excludeperiods', excludeperiods,'minstd',minthresh);
end

%Go through each cell and determine if it is phase locked to ripples and
%gammas during ripple events and out of ripple events
if ~isempty(riptimes)
    %Exclude ripples occuring close together for preceding modulationg
    valid_ripples = [1000; riptimes(2:end,2)-riptimes(1:end-1,1)];
    valid_ripples = valid_ripples > 1;

    for cellcount = 1:size(index,1)
        ind = index(cellcount,:);

        if ~isempty(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data)
            spiketimes = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);
            
            if isfield(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)},'lowgammaphase')
                gam = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.lowgammaphase;
            else
                gam = nan(size(spiketimes));
            end
            if isfield(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)},'ripplephase')
                rip = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.ripplephase;
            else
                rip = nan(size(spiketimes));
            end
        else
            spiketimes = [];;
            gam = [];
            rip = [];
        end

        %Find valid spikes
        spikebins = periodAssign(spiketimes, riptimes(valid_ripples,[1 2]));
        maxthresh = zeros(size(spikebins));
        maxthresh(spikebins>0) = ripstd(spikebins(spikebins>0));
        
        if ~isempty(spiketimes)
            validspikes = find(~isExcluded(spiketimes, excludeperiods) & spikebins);
            spikebins = spikebins(validspikes);
            spiketimes = spiketimes(validspikes);
            maxthresh = maxthresh(validspikes);
            
            gam = gam(validspikes);
            rip = rip(validspikes);
            tmpcelldata = [spiketimes spikebins maxthresh gam rip];
            tmpcelldata(:,6) = cellcount;
           
        else
            tmpcelldata = [0 0 NaN NaN NaN cellcount];
        end

        %add cell location 1 if CA1, 3 if CA3
        if ~isfield(cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)},'area')
            tmpcelldata(:,7) = 0;
        elseif isequal('CA1',cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)}.area)
            tmpcelldata(:,7) = 1;
        elseif isequal('CA3',cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)}.area)
            tmpcelldata(:,7) = 3;
        else
            tmpcelldata(:,7) = 0;
        end

        celldata = [celldata; tmpcelldata];

        spikecount = zeros(1,size(riptimes(valid_ripples),1));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        spikecounts = [spikecounts; spikecount];
    end

    
    %Sort all spikes by time
    celldata = sortrows(celldata,1);

    celldata(:,8) = 0;
    valid = find(sum(spikecounts>0)>=min_cells);
    for i = 1:length(valid)
        celldata(celldata(:,2)==valid(i),8) = 1;
    end
    out.celldata = celldata(logical(celldata(:,8)),[1 2 3 4 5 6 7]);
else
    out.celldata = [];
end

end
