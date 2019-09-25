function out = calcspikingmodinripples(index, excludeperiods, spikes, ripples,cellinfo, varargin)
% Computes the modulation of spikes by gamma during ripples.
% out.celldata is a Nx6 matrix with N rows for each spike and 6 rows
%   1: spiketime
%   2: which ripple the spike occured in
%   3: gamma phase on the local tetrode
%   4: gamma phase on the reference tetrode
%   5: which cell
%   6: area recorded

% out.celldata_first is a Nx6 matrix with N rows for each spike and 6 rows
%   corresponding to the information for the first spike
% out.preceding is an Nx6 matrix with N rows for all spikes that occured in
%   the 500ms before each ripple

%Set options
cellfilter = [];
min_cells = 1;

celldata = []; celldata_first =[]; preceding = []; phase_all = nan(size(index,1),2);
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
            case 'frequency'
                freq = varargin{option+1};
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
    riptimes = getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
        '(isequal($area, ''CA1''))',...
        'excludeperiods', excludeperiods,'minstd',minthresh);
else
    riptimes = getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
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
            if isfield(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)},'globalgammaphase')
                rgam = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.globalgammaphase;
            else
                rgam = nan(size(spiketimes));
            end
        else
            spiketimes = [];
            gam = [];
            rgam = [];
        end

        %Determine preferred phase for all time
        if ~isempty(rgam)
            tmp = rayleigh_test(rgam);
            phase_all(cellcount,1) = tmp.theta;
        end
        %Find valid spikes
        spikebins = periodAssign(spiketimes, riptimes(valid_ripples,[1 2]));
        precedingbins = periodAssign(spiketimes, [riptimes(valid_ripples,1)-0.5 riptimes(valid_ripples,1)]);
        
        pspiketimes = spiketimes;   pgam = gam; prgam = rgam;
        
        if ~isempty(spiketimes)
            validspikes = find(~isExcluded(spiketimes, excludeperiods) & spikebins);
            spikebins = spikebins(validspikes);
            spiketimes = spiketimes(validspikes);

            gam = gam(validspikes);
            rgam = rgam(validspikes);
            tmpcelldata = [spiketimes spikebins gam rgam];
            tmpcelldata(:,5) = cellcount;
            
            validspikes = find(~isExcluded(pspiketimes, excludeperiods) & precedingbins);
            precedingbins = precedingbins(validspikes);
            pspiketimes = pspiketimes(validspikes);
            pgam = pgam(validspikes);
            prgam = prgam(validspikes);
            tmppreceding = [pspiketimes precedingbins pgam prgam];
            tmppreceding(:,5) = cellcount;
           
            %Identify first spike of each cell
            if ~isempty(spiketimes)
                valid = logical([1; diff(spikebins)>0]);
                spikebins_first = spikebins(valid);
                spiketimes_first = spiketimes(valid);
                gam_first = gam(valid);
                rgam_first = rgam(valid);
                
                tmpcelldata_first = [spiketimes_first spikebins_first gam_first rgam_first];
                tmpcelldata_first(:,5) = cellcount;
            else
                tmpcelldata_first = [0 0 NaN NaN cellcount];
            end
            
        else
            tmpcelldata = [0 0 NaN NaN cellcount];
            tmppreceding = [0 0 NaN NaN cellcount];
        end

        %add cell location 1 if CA1, 3 if CA3
        if ~isfield(cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)},'area')
            tmpcelldata(:,6) = 0;
            tmppreceding(:,6) = 0;
            tmpcelldata_first(:,6) = 0;
        elseif isequal('CA1',cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)}.area)
            tmpcelldata(:,6) = 1;
            tmppreceding(:,6) = 1;
            tmpcelldata_first(:,6) = 1;
            phase_all(cellcount,2) = 1;
        elseif isequal('CA3',cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)}.area)
            tmpcelldata(:,6) = 3;
            tmppreceding(:,6) = 3;
            tmpcelldata_first(:,6) = 3;
            phase_all(cellcount,2) = 3;
        else
            tmpcelldata(:,6) = 0;
            tmpcelldata_first(:,6) = 0;
        end

        celldata = [celldata; tmpcelldata];
        preceding = [preceding; tmppreceding];
        
        celldata_first = [celldata_first; tmpcelldata_first];
        spikecount = zeros(1,size(riptimes,1));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        spikecounts = [spikecounts; spikecount];
    end

    %Sort all spikes by time
    celldata = sortrows(celldata,1); celldata_first = sortrows(celldata_first,1);
    preceding = sortrows(preceding,1);

    celldata(:,7) = 0; celldata_first(:,7) = 0;
    valid = find(sum(spikecounts>0)>=min_cells);
    for i = 1:length(valid)
        celldata(celldata(:,2)==valid(i),7) = 1;
        celldata_first(celldata_first(:,2)==valid(i),7) = 1;
    end
    out.celldata = celldata(logical(celldata(:,7)),[1 2 3 4 5 6]);
    out.preceding = preceding;
    out.celldata_first = celldata_first(logical(celldata_first(:,7)),[1 2 3 4 5 6]);
    out.phase_all = phase_all;
else
    out.celldata = [];
    out.preceding = [];
    out.celldata_first = [];
    out.phase_all = [];
end

end
