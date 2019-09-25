function out = calcpopulationlinfields(indices, excludeperiods, spikes, linpos, binsize, peakthresh)
%trajdata = FILTERCALCLINFIELDS(index, excludeperiods, spikes, linpos)
%trajdata = FILTERCALCLINFIELDS(index, excludeperiods, spikes, linpos, binsize)
%
%Calculates the linear occupancy normalized firing rate for all the cells.
%
%spikes - the 'spikes' cell array for the day you are analyzing
%linpos - the output of LINEARDAYPROCESS for the day you are analyzing. 
%indices - [day epoch tetrode cell]
%binsize- the length of each spatial bin (default 2cm)
%excludeperiods - [start end] times for each exlcude period
%
%The output is a structure. 


warning('OFF','MATLAB:divideByZero');
if (nargin < 5)
    binsize = 2;
end

%Define linear position information
index = [indices(1,1) indices(1,2)];
statematrix = linpos{index(1)}{index(2)}.statematrix;
wellSegmentInfo = linpos{index(1)}{index(2)}.wellSegmentInfo;
segmentInfo = linpos{index(1)}{index(2)}.segmentInfo;
trajwells = linpos{index(1)}{index(2)}.trajwells;
statevector = statematrix.traj;
lindist = statematrix.lindist;
lindist = gaussSmooth(lindist, 30);

%Apply excludeperiods to find valid trajectories
if iscell(statevector)
    statevector = statevector{6};
end
statevector(find(isExcluded(statematrix.time, excludeperiods))) = -1;

%calculate information about the track and the behavior
for i = 1:size(trajwells,1)
    %calculate the linear length of the trajectory
    trajlength(i) = sum(segmentInfo.segmentLength(wellSegmentInfo.pathTable{trajwells(i,1),wellSegmentInfo.segmentIndex(trajwells(i,2))}));
end

try
    %Define the time, distance, and trajectory for all good times.
    goodlocationind = (find(statevector ~= -1 ));
    goodlocations = [statematrix.time(goodlocationind) lindist(goodlocationind) statevector(goodlocationind)];
end

% Initialize output
out.traj = [];
out.dist = [];
out.homesegmentlength = [];
out.index = [];
out.rates = [];
out.posprob = [];

%go through each cell and calculate the linearized rates
for cellcount = 1:size(indices,1)
    
    index = indices(cellcount,:);
    cellspikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;

    trajnum = max(statevector);
    timestep = statematrix.time(2,1) - statematrix.time(1,1);

    posindexfield = 7;
    if ~isempty(cellspikes)
        cellspikes = cellspikes(:,[1 posindexfield]);
        %Add the trajectory for each spike
        cellspikes(:,3) = statevector(cellspikes(:,2));
        %Add the linear distance for each spike
        cellspikes(:,4) = lindist(cellspikes(:,2));
    else
        cellspikes = [0 0 -1 0];
    end

    goodspikes = [];
    goodspikeind = (cellspikes(:,3) ~= -1);
    goodspikes = cellspikes(goodspikeind,:);
   
    tmprates = [];
    tmptraj = [];
    tmprates = [];
    tmpdist = [];
    tmpocc = [];
    for i = 1:trajnum

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
                %the 2nd column of binvector is occupancy
                binvector(j,2) = sum(binind == j) * timestep;
                %the 3rd column of binvector is the spike count
                binvector(j,3) = sum(spikebinind == j);
            end
            %the 4th column of binvector is occ-normalized firing rate
            binvector(:,4) = binvector(:,3)./binvector(:,2);
            nonfinite = find(~isfinite(binvector(:,4)));
            %the 6th column of binvector is the smoothed occupancy
            binvector(:,6) = gaussSmooth(binvector(:,2),2);
            %the 7th column of binvector is the smoothed spike count
            binvector(:,7) = gaussSmooth(binvector(:,3),2);
            %the 5th column of binvector is the smoothed occ-normalized
            %firing rate (with the baseline rate added)
            binvector(:,5) = (binvector(:,7)./binvector(:,6))+.05;
            %Exclued bins with very low occupancy and replace with NaN
            lowocc = find(binvector(:,6) < binsize*.1);
            binvector(nonfinite,4) = nan;
            binvector(lowocc,5) = nan;

            tmprates = [tmprates binvector(:,5)'];
            tmptraj = [tmptraj ones(1,size(binvector,1))*i];
            tmpdist = [tmpdist binvector(:,1)'];
            tmpocc = [tmpocc binvector(:,6)'];
            %trajdata{i} = binvector;
        end
    end

    if max(tmprates >= peakthresh)
        out.dist = tmpdist;
        out.traj = tmptraj;
        out.homesegmentlength = segmentInfo.segmentLength(1);
% WHAT IS THE GOAL OF THE POSPROB FIELD?
        out.posprob = tmpocc/noNanSum(tmpocc);
        out.rates = [out.rates;tmprates];
        out.index = [out.index; index];
    end
end

warning('ON','MATLAB:divideByZero');



function  out = gaussSmooth(vector, binrange)

paddinglength = round(binrange*2.5);
padding = ones(paddinglength,1);

out = smoothvect([padding*vector(1) ; vector; padding*vector(end)],gaussian(binrange,binrange*7));
out = out(paddinglength+1:end-paddinglength);



