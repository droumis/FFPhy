function out = calcswrrate(index, excludeperiods, spikes, ripples,cellinfo, varargin)
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

celldata = []; preceding = [];
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

% Separate cells that ever have a place field with non-place fields
non_place_fields = evaluatefilter(cellinfo{index(1,1)},'$peakrate>3');
non_place_fields = rowfind(index(:,3:4),non_place_fields(:,2:end))==0;

index(:,5) = 0;
index(non_place_fields,5) = 1;

%Go through each cell and determine if it is phase locked to ripples and
%gammas during ripple events and out of ripple events
if ~isempty(riptimes)
    %Exclude ripples occuring close together for preceding modulationg
    valid_ripples = [1000; riptimes(2:end,2)-riptimes(1:end-1,1)];
    valid_ripples = valid_ripples > 1;
    
    %Find total time
    preceding_time = length(riptimes(valid_ripples,2))*0.5;

    for cellcount = 1:size(index,1)
        ind = index(cellcount,:);

        if ~isempty(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data)
            spiketimes = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);
        else
            spiketimes = [];
        end
        pspiketimes = spiketimes;
        
        %Find valid spikes
        spikebins = periodAssign(spiketimes, riptimes(:,[1 2]));
        precedingbins = periodAssign(spiketimes, [riptimes(valid_ripples,1)-0.5 riptimes(valid_ripples,1)]);
    
        if ~isempty(spiketimes)
            validspikes = find(~isExcluded(spiketimes, excludeperiods) & spikebins);
            spikebins = spikebins(validspikes);
            spiketimes = spiketimes(validspikes);
            tmpcelldata = [spiketimes spikebins];
            tmpcelldata(:,3) = cellcount;
            tmpcelldata(:,4) = ind(5);
            validspikes = find(~isExcluded(pspiketimes, excludeperiods) & precedingbins);
            precedingbins = precedingbins(validspikes);
            pspiketimes = pspiketimes(validspikes);
            tmppreceding = [pspiketimes precedingbins];
            tmppreceding(:,3) = cellcount;
            tmppreceding(:,4) = ind(5);
           
        else
            tmpcelldata = [0 0 cellcount ind(5)];
            tmppreceding = [0 0 cellcount ind(5)];
            
        end

        %add cell location 1 if CA1, 3 if CA3
        if ~isfield(cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)},'area')
            tmpcelldata(:,5) = 0;
            tmppreceding(:,5) = 0;
        elseif isequal('CA1',cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)}.area)
            tmpcelldata(:,5) = 1;
            tmppreceding(:,5) = 1;
        elseif isequal('CA3',cellinfo{ind(1)}{ind(2)}{ind(3)}{ind(4)}.area)
            tmpcelldata(:,5) = 3;
            tmppreceding(:,5) = 3;
        else
            tmpcelldata(:,5) = 0;
            tmpcelldata_first(:,5) = 0;
        end

        celldata = [celldata; tmpcelldata];
        preceding = [preceding; tmppreceding];
        
        spikecount = zeros(1,size(riptimes,1));
        for i = 1:length(spikebins)
            spikecount(spikebins(i)) = spikecount(spikebins(i))+1;
        end
        spikecounts = [spikecounts; spikecount];
    end

    %Calc mean rate for preceding time
    place = preceding(:,4)==0 & preceding(:,5) == 1;
    ca1_preceding_place_cells = sum(place)./preceding_time;
    ca1_preceding_place_cells_nspikes = sum(place);
    
    place = preceding(:,4)==0 & preceding(:,5)==3;
    ca3_preceding_place_cells = sum(place)./preceding_time;
    ca3_preceding_place_cells_nspikes = sum(place);
    
    place = preceding(:,4)==1 & preceding(:,5) == 1;
    ca1_preceding_not_place_cells = sum(place)./preceding_time;
    ca1_preceding_not_place_cells_nspikes = sum(place);
    
    place = preceding(:,4)==1 & preceding(:,5)==3;
    ca3_preceding_not_place_cells = sum(place)./preceding_time;
    ca3_preceding_not_place_cells_nspikes = sum(place);

    
    %Calc mean SWR rate
    celldata = sortrows(celldata,2);
    celldata(:,6) = 0;
    valid = find(sum(spikecounts>0)>=min_cells);
    for i = 1:length(valid)
        celldata(celldata(:,2)==valid(i),6) = 1;
    end
    swr_time = sum(riptimes(valid,2)-riptimes(valid,1));
    
    %Calc mean rate for SWRs
    place = celldata(:,4)==0 & celldata(:,5) == 1 & celldata(:,6)==1;
    ca1_place_cells = sum(place)./swr_time;
    ca1_place_cells_nspikes = sum(place);
    
    place = celldata(:,4)==0 & celldata(:,5) == 3 & celldata(:,6)==1;
    ca3_place_cells = sum(place)./swr_time;
    ca3_place_cells_nspikes = sum(place);
    
    place = celldata(:,4)==1 & celldata(:,5) == 1 & celldata(:,6)==1;
    ca1_not_place_cells = sum(place)./swr_time;
    ca1_not_place_cells_nspikes = sum(place);
    
    place = celldata(:,4)==1 & celldata(:,5) == 3 & celldata(:,6)==1;
    ca3_not_place_cells = sum(place)./swr_time;
    ca3_not_place_cells_nspikes = sum(place);
    
    ca1 = length(unique(celldata(celldata(:,5)==1 & celldata(:,6)==1,3)));
    ca3 = length(unique(celldata(celldata(:,5)==3 & celldata(:,6)==1,3)));
    
    out = [ca1 ca1_place_cells ca1_place_cells_nspikes ca1_not_place_cells ca1_not_place_cells_nspikes ...
        ca1_preceding_place_cells ca1_preceding_place_cells_nspikes ca1_preceding_not_place_cells ca1_preceding_not_place_cells_nspikes ...
        ca3 ca3_place_cells ca3_place_cells_nspikes ca3_not_place_cells ca3_not_place_cells_nspikes ...
        ca3_preceding_place_cells ca3_preceding_place_cells_nspikes ca3_preceding_not_place_cells ca3_preceding_not_place_cells_nspikes];

else
    out = [];
end

end