function [output] = openfieldratemap(spikes,pos, binsize, std)

% See also calcopenfieldoccupancy and twodoccupancy
%[output] = openfieldratemap(spikes,pos, binsize, std)
%[output] = openfieldratemap(spikes,pos, binsize)
%
%Calculates the 2d occupancy normalized firing rate for the cell
%
%spikes - the spike data structure 
%pos - the position data structure
%binsize- the length of each spatial bin (default 1 cm)
%std - defines the shape of the 2d gaussian used to smooth spikerate.
%              (default 1)
%
%The output is a structure with n matrices. The matrices are: occupancy, 
%bin vector x, bin vector y, bin spike count, occ normailized firing per bin, and smoothed
%occ normalized firing. 
%

warning('OFF','MATLAB:divideByZero');


if (nargin < 4)
    std = 1;
end

if (nargin < 3)
    binsize = 1;
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
    minx = floor(min(tmpposition(:,1))) - 1;
    maxx = ceil(max(tmpposition(:,1))) + 1;
    binx = (minx:binsize:maxx);
    miny = floor(min(tmpposition(:,2))) - 1;
    maxy = ceil(max(tmpposition(:,2))) + 1;
    biny = (miny:binsize:maxy);
    %[output.occupancy output.xticks output.yticks] = hist2(tmpposition, binx, biny);
    [output.occupancy output.xticks output.yticks] = HIST2(tmpposition(:,1), tmpposition(:,2), binx, biny);
    %[output.occupancy output.xticks output.yticks] = hist2d(tmpposition, min([length(binx),length(biny)]));
    
    nonzero = find(output.occupancy ~= 0);
    %[output.spikes BX BY] = hist2(tmpspikes, binx, biny);
    [output.spikes BX BY] = HIST2(tmpspikes(:,1), tmpspikes(:,2), binx, biny);
    %[output.spikes BX BY] = hist2d(tmpspikes, binsize);
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

warning('ON','MATLAB:divideByZero');
