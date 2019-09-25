function [output] = twodoccupancy(spikes,statevector, linpos, pos, index, binsize, std)
%[output] = twodoccupancy(spikes,statevector, linpos, pos, index, binsize)
%[output] = twodoccupancy(spikes,statevector, linpos, pos, index)
%
%Calculates the 2d occupancy normalized firing rate for the cell and
%organizes the output into the different trajectories.
%
%spikes - the 'spikes' cell array for the day you are analyzing
%linpos - the output of LINEARDAYPROCESS for the day you are analyzing.
%pos - the output of nspike_fixpos
%index - [day epoch tetrode cell]
%statevector - the outputs of GETBEHAVESTATE. This is a vector with the traj
%              number for each position time (1 based). -1 values signify
%              invalid times and are not used.
%binsize- the length of each spatial bin (default 1 cm)
%std - defines the shape of the 2d gaussian used to smooth spikerate.
%              (default 1)
%
%The output is a cell array where each cell contains a cell
%descibing one trajectory.  These cells contain n matrices. The matrices are: occupancy, 
%bin vector x, bin vector y, bin spike count, occ normailized firing per bin, and smoothed
%occ normalized firing. If the cell is empty, the animal did not enter that trajectory.
%

warning('OFF','MATLAB:divideByZero');

statematrix = linpos{index(1)}{index(2)}.statematrix;
spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
spikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
pos = pos{index(1)}{index(2)}.data;

if (nargin < 7)
    std = 1;
else
    %std = user defined;
end

if (nargin < 6)
    binsize = 1;
else
    %binsize = user defined;
end

posx = 2;
posy = 3;
posindexfield = 7;

if ~isempty(spikes)
    spikes = spikes(:,[1 posx posy posindexfield]);
    spikes(:,5) = statevector(spikes(:,4)); %add the traj state for each spike
    %posindexfield = isdatafield(spikesfields,'posindex');

else
    spikes = [0 0 -1];
end

trajnum = max(statevector);

    timestep = statematrix.time(2,1) - statematrix.time(1,1);
    
    goodspikes = [];

    goodspikeind = (spikes(:,5) ~= -1);
    %create a list of the non sharp-wave spikes [time traj linearloc]
    goodspikes = spikes(goodspikeind,:);
    %make a cell array, where each cell contains data for one trajectory.
    %inside each cell [binlocation occupancy spikecount firingrate]
    
    try 
        trajdata = cell(1,trajnum);
        goodlocationind = (find(statevector ~= -1 ));
        goodlocations = [statematrix.time(goodlocationind) pos(goodlocationind,[2 3]) statevector(goodlocationind)]; %CHANGED the positionindex at all valid times
    end

%     if isempty(goodspikes) %there were no useable spikes
%         goodspikes = [];
%         trajdata = [];
%     end

for i = 1:length(trajdata)
   
    %get all the position indices and spikes when the animal was on the ith
    %trajectory

    tmpposition = goodlocations(find(goodlocations(:,4) == i),[2 3]);
    tmpspikes = goodspikes(find(goodspikes(:,5) == i),[2 3]);
    
    if ~isempty(tmpposition)
        minx = floor(min(tmpposition(:,1)));
        maxx = ceil(max(tmpposition(:,1)));
        binx = (minx:binsize:maxx);
        miny = floor(min(tmpposition(:,2)));
        maxy = ceil(max(tmpposition(:,2)));
        biny = (miny:binsize:maxy);
        [output(i).occupancy output(i).xticks output(i).yticks] = hist2(tmpposition, binx, biny);
        
        if ~isempty(tmpspikes)
        [output(i).spikes BX BY] = hist2(tmpspikes, binx, biny);
        nonzero = find(output(i).occupancy ~= 0);
        output(i).spikerate = zeros(size(output(i).spikes));
        output(i).spikerate(nonzero) = output(i).spikes(nonzero) ./(timestep* output(i).occupancy(nonzero) );

%smoothed occupancy
        g = gaussian2(std,(6*std));
        output(i).smoothedspikerate = filter2(g,(output(i).spikerate));
        smoothedoccupancy = []
        smoothedoccupancy = zeros(size(output(i).spikes));
        smoothedoccupancy = filter2(g, output(i).occupancy);
        zero = find(smoothedoccupancy == 0);
        output(i).smoothedspikerate(zero) = -1;
        end
    end
end

warning('ON','MATLAB:divideByZero');