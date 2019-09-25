function out = DFAsj_filtercalclinfields_tf(index, excludeperiods, spikes, linpos, varargin)

% Equivalent to sj_calclinfields and DFsj_filtercalclinfields - WITH TIMEFILTER
% DO NOT NEED POS for this version UNLIKE DFsj_filtercalclinfields since
% speed is supposed to have been filtered out by timefilter
%
%trajdata = FILTERCALCLINFIELDS(index, excludeperiods, spikes, linpos, pos)
%trajdata = FILTERCALCLINFIELDS(index, excludeperiods, spikes, linpos, pos, options)
%
%Calculates the linear occupancy normalized firing rate for the cell and
%organizes the output into the different trajectories.
%
%spikes - the 'spikes' cell array for the day you are analyzing
%linpos - the output of LINEARDAYPROCESS for the day you are analyzing. 
%index - [day epoch tetrode cell]
%binsize- the length of each spatial bin (default 2cm)
%excludeperiods - [start end] times for each exlcude period
%
% options are
%	'binsize'	binsize (default 2 cm)
%	'phasedist'	0 or 1 	Adds a phasedist field for each trajectory
%				that contains the phase and linear distance for 
%				each spike on that trajectory (default 0)
%
%The output is a structure where the trajectory field contains a matrices
%descibing the trajectories.  These matrices are n by 7, where n is the
%number of spatial bins. The columns are: linear bin location, bin
%occupancy (seconds), bin spike count, occ normailized firing per bin, and
%smoothed occ normalized firing, smoothed occupancy, smoothed spike count. If the cell is empty, the animal did not
%enter that trajectory.
%

%day = index(1)

% parse the options
phasedist = 0;
binsize = 2;


for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'binsize'
                binsize = varargin{option+1};
            case 'phasedist'
                phasedist = varargin{option+1};
	    otherwise
                error(['Option ',varargin{option},' unknown.']);
	end   
    else
        error('Options must be strings, followed by the variable');
    end
end

index;
warning('OFF','MATLAB:divideByZero');
statematrix = linpos{index(1)}{index(2)}.statematrix;
wellSegmentInfo = linpos{index(1)}{index(2)}.wellSegmentInfo;
segmentInfo = linpos{index(1)}{index(2)}.segmentInfo;
trajwells = linpos{index(1)}{index(2)}.trajwells;
spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
spikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;

lindist = statematrix.lindist;
statevector = statematrix.traj;


% Use Exclude Periods for TimeFilter version in addition to statevector=-1
% that already exists based on behavestate
statevector(find(isExcluded(statematrix.time, excludeperiods))) = -1; % Based on exclude time, rather than inputing statevector as in calclinfields


if (length(segmentInfo.segmentLength) > 1)
    for i = 1:size(trajwells,1)
        %calculate the linear length of the trajectory
        trajlength(i) = sum(segmentInfo.segmentLength(wellSegmentInfo.pathTable{trajwells(i,1),wellSegmentInfo.segmentIndex(trajwells(i,2))}));
    end
else
    trajlength(1) = segmentInfo.segmentLength(1);
end


%posindexfield = isdatafield(spikesfields,'posindex');
posindexfield = 7;
if ~isempty(spikes)
    tmpspikes = spikes;
    spikes = spikes(:,[1 posindexfield]);
    % Find any posidxs that are outside the range of linpos (statevector)
    outliers = find(spikes(:,2) > length(statevector));
    spikes(outliers,:)=[];
    spikes(:,3) = statevector(spikes(:,2)); %add the traj state for each spike
    spikes(:,4) = lindist(spikes(:,2)); %add the linear distance
    if (phasedist)
        spikes(:,5) = tmpspikes(:,5); % theta phase
    end
else
    spikes = [0 0 -1 0 0];
end


trajnum = max(statevector);
timestep = statematrix.time(2,1) - statematrix.time(1,1);
goodspikes = [];
goodspikeind = (spikes(:,3) ~= -1);

%create a list of the non sharp-wave spikes [time traj linearloc]
goodspikes = spikes(goodspikeind,:);
%make a cell array, where each cell contains data for one trajectory.
%inside each cell [binlocation occupancy spikecount firingrate]

trajdata = cell(1,trajnum);
goodlocationind = find((statevector ~= -1));
goodlocations = [statematrix.time(goodlocationind) lindist(goodlocationind) statevector(goodlocationind)]; %the linear locations at all valid times


trajval = [];
for i = 1:length(trajdata)

    %get all the linear locations when the animal was on the ith
    %trajectory
    tmplinloc = goodlocations(find(goodlocations(:,3) == i),2);  % lindist for current traj
    tmpspikes = goodspikes(find(goodspikes(:,3) == i),:); % Spikes for current traj
    if ~isempty(tmplinloc)
        findtrajnum = (i+rem(i,2))/2;
        minloc = 0;
        maxloc = trajlength(findtrajnum);

        
        tmpvec = [minloc:binsize:maxloc];
        binvector = zeros(length(tmpvec),5);
        %the 1st column of binvector is the location for that bin
        binvector(:,1) = tmpvec(1:end);
        %find which bins the linear locations and spikes fall into
        binind = lookup(tmplinloc,binvector(:,1));
        if ~isempty(tmpspikes)
            spikebinind = lookup(tmpspikes(:,4),binvector(:,1));
        else
            spikebinind = [];
        end
        %sum up the occupancy and spikes in each bin
        for j = 1:size(binvector,1)
            binvector(j,2) = sum(binind == j) * timestep;  %the second column is occupancy
            binvector(j,3) = sum(spikebinind == j); %3rd column is spike count
        end
        binvector(:,4) = binvector(:,3)./binvector(:,2); %4th column is occ-normalized firing rate
        nonfinite = find(~isfinite(binvector(:,4)));
        binvector(:,6) = gaussSmooth(binvector(:,2),2); %6th column is smoothed occupancy
        binvector(:,7) = gaussSmooth(binvector(:,3),2); %7th column is smoothed spike count
        binvector(:,5) = binvector(:,7)./binvector(:,6); %5th column is smoothed occ-normalized firing rate
        lowocc = find(binvector(:,6) < binsize*.1); %very low occupancy bins should be excluded


        %after the smoothing, turn the firing rate bins that had low occupancy to nan's.
        binvector(nonfinite,4) = nan;
        binvector(lowocc,5) = nan;

        trajdata{i} = binvector;
	if (phasedist) 
	    out.phasedist{i}.phase = tmpspikes(:,5);
	    out.phasedist{i}.dist = tmpspikes(:,4);
	end
    end

end

out.trajdata = trajdata;
out.index = index;
warning('ON','MATLAB:divideByZero');



function  out = gaussSmooth(vector, binrange)

paddinglength = round(binrange*2.5);
padding = ones(paddinglength,1);

out = smoothvect([padding*vector(1) ; vector; padding*vector(end)],gaussian(binrange,binrange*5));
out = out(paddinglength+1:end-paddinglength);



