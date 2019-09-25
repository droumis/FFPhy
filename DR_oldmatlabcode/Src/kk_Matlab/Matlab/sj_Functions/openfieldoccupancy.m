function [output] = openfieldoccupancy(spikes,pos, index, binsize, std)
%[output] = openfieldoccupancy(spikes,pos, index, binsize, std)
%[output] = openfieldoccupancy(spikes,pos, index, binsize)
%
%Calculates the 2d occupancy normalized firing rate for the cell and
%organizes the output into the different trajectories.
%
%spikes - the 'spikes' cell array for the day you are analyzing
%pos - the output of nspike_fixpos
%index - [day epoch tetrode cell]
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

spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
spikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
pos = pos{index(1)}{index(2)}.data;

if (nargin < 5)
    std = 1;
else
    %std = user defined;
end

if (nargin < 4)
    binsize = 1;
else
    %binsize = user defined;
end

posx = 2;
posy = 3;
posindexfield = 7;

if ~isempty(spikes)
    spikes = spikes(:,[1 posx posy posindexfield]);

else
    spikes = [0 0 -1];
end

    timestep = pos(2,1) - pos(1,1);
    
    goodspikes = [];
    goodspikes = spikes;
    %make a cell array, where each cell contains data for one trajectory.
    %inside each cell [binlocation occupancy spikecount firingrate]

    tmpposition = (pos(:,[2 3]));
    tmpspikes = (goodspikes(:,[2 3]));
    
    if ~isempty(tmpposition)
        minx = floor(min(tmpposition(:,1)));
        maxx = ceil(max(tmpposition(:,1)));
        binx = (minx:binsize:maxx);
        miny = floor(min(tmpposition(:,2)));
        maxy = ceil(max(tmpposition(:,2)));
        biny = (miny:binsize:maxy);
        [output.occupancy output.xticks output.yticks] = hist2(tmpposition, binx, biny);
        
        if ~isempty(tmpspikes)
        [output.spikes BX BY] = hist2(tmpspikes, binx, biny);
        nonzero = find(output.occupancy ~= 0);
        output.spikerate = zeros(size(output.spikes));
        output.spikerate(nonzero) = output.spikes(nonzero) ./(timestep* output.occupancy(nonzero) );

%smoothed occupancy
        g = gaussian2(std,(6*std));
        output.smoothedspikerate = filter2(g,(output.spikerate)); % is this the right filter?
        smoothedoccupancy = [];
        smoothedoccupancy = zeros(size(output.spikes));
        smoothedoccupancy = filter2(g, output.occupancy);
        zero = find(smoothedoccupancy == 0);
        output.smoothedspikerate(zero) = -1;
        end
    end

warning('ON','MATLAB:divideByZero');