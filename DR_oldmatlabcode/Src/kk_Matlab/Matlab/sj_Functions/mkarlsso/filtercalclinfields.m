function out = filtercalclinfields(index, excludeperiods, spikes, linpos, binsize)
%trajdata = FILTERCALCLINFIELDS(index, excludeperiods, spikes, linpos)
%trajdata = FILTERCALCLINFIELDS(index, excludeperiods, spikes, linpos, binsize)
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
%The output is a structure where the trajectory field contains a matrices
%descibing the trajectories.  These matrices are n by 5, where n is the
%number of spatial bins. The columns are: linear bin location, bin
%occupancy (seconds), bin spike count, occ normailized firing per bin, and
%smoothed occ normalized firing. If the cell is empty, the animal did not
%enter that trajectory.
%



warning('OFF','MATLAB:divideByZero');
statematrix = linpos{index(1)}{index(2)}.statematrix;
wellSegmentInfo = linpos{index(1)}{index(2)}.wellSegmentInfo;
segmentInfo = linpos{index(1)}{index(2)}.segmentInfo;
trajwells = linpos{index(1)}{index(2)}.trajwells;
spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
spikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
statevector = statematrix.traj;
lindist = statematrix.lindist;
statevector(find(isExcluded(statematrix.time, excludeperiods))) = -1;

for i = 1:size(trajwells,1)
    %calculate the linear length of the trajectory
    trajlength(i) = sum(segmentInfo.segmentLength(wellSegmentInfo.pathTable{trajwells(i,1),wellSegmentInfo.segmentIndex(trajwells(i,2))}));
end
if (nargin < 5)
    binsize = 2;
end

%posindexfield = isdatafield(spikesfields,'posindex');
posindexfield = 7;
if ~isempty(spikes)
    spikes = spikes(:,[1 posindexfield]);
    spikes(:,3) = statevector(spikes(:,2)); %add the traj state for each spike
    spikes(:,4) = lindist(spikes(:,2)); %add the linear distance

else
    spikes = [0 0 -1 0];
end

trajnum = max(statevector);
timestep = statematrix.time(2,1) - statematrix.time(1,1);
goodspikes = [];

goodspikeind = (spikes(:,3) ~= -1);
%create a list of the non sharp-wave spikes [time traj linearloc]
goodspikes = spikes(goodspikeind,:);
%make a cell array, where each cell contains data for one trajectory.
%inside each cell [binlocation occupancy spikecount firingrate]
try
    trajdata = cell(1,trajnum);
    goodlocationind = (find(statevector ~= -1 ));
    goodlocations = [statematrix.time(goodlocationind) lindist(goodlocationind) statevector(goodlocationind)]; %the linear locations at all valid times
end

trajval = [];
for i = 1:length(trajdata)

    %get all the linear locations when the animal was on the ith
    %trajectory
    tmplinloc = goodlocations(find(goodlocations(:,3) == i),2);
    tmpspikes = goodspikes(find(goodspikes(:,3) == i),:);
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


