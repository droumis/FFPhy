function [out] = calcopenfieldoccupancy(index, excludetimes, spikes,pos, varargin)
%[out] = openfieldoccupancy(index, excludetimes, spikes,pos, options)
%
%Calculates the 2d occupancy normalized firing rate for the cell.  Does NOT
%separate by trajectory (to do so, see twodoccupancy)
%
% excludetime -- times want to exclude from analysis
% spikes - the 'spikes' cell array for the day you are analyzing
% pos - the output of nspike_fixpos
% index - [day epoch tetrode cell]
% options:
%  binsize- the length of each spatial bin (default 1 cm)
%  std - defines the shape of the 2d gaussian used to smooth spikerate.
%              (default 1)
%  appendindex - 0 or 1, 1 appends index infront of output (default 0 )
%
%The output is a structure with fields: occupancy, bin vector x (.xticks), bin vector
%y (.yticks), bin spike count (.spikes), occ normailized firing per bin (.spikerate), and smoothed occ
% normalized firing (.smoothedspikerate).


appendindex = 0;
std = 1;
binsize = 1;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'std'
                std = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end
warning('OFF','MATLAB:divideByZero');

spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
spikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
pos = pos{index(1)}{index(2)}.data;

if (nargin < 6)
    std = 1;
else
    %std = user defined;
end

if (nargin < 5)
    binsize = 1;
else
    %binsize = user defined;
end

posx = 2;
posy = 3;
posindexfield = 7;

if ~isempty(spikes)
    spikes = spikes(:,[1 posx posy posindexfield]); %columns: time, x y, posindex

else
    spikes = [0 0 -1];
end

timestep = mean(diff(pos(:,1)));;
%filter out exclude times
indgoodspikes = ~isExcluded(spikes(:,1), excludetimes);
goodspikes = spikes(indgoodspikes,:);  %select spikes not excluded by exclude times
indgoodpos =  ~isExcluded(pos(:,1), excludetimes);
goodpos = pos(indgoodpos,:);  %select spikes not excluded by exclude times

tmpposition = (goodpos(:,[2 3])); %xy position
tmpspikes = (goodspikes(:,[2 3]));

if ~isempty(tmpposition)
    minx = floor(min(tmpposition(:,1)));
    maxx = ceil(max(tmpposition(:,1)));
    binx = (minx:binsize:maxx);
    miny = floor(min(tmpposition(:,2)));
    maxy = ceil(max(tmpposition(:,2)));
    biny = (miny:binsize:maxy);
    [output.occupancy output.xticks output.yticks] = hist2(tmpposition(:,1), tmpposition(:,2), binx, biny);

    if ~isempty(tmpspikes)
        [output.spikes BX BY] = hist2(tmpspikes(:,1), tmpspikes(:,2), binx, biny);
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
        output.smoothedspikerate(zero) = -2; %no occupancy is negative and diff/darker from occupancy but no spikes
    end
end

if appendindex == 0
    out = output;
elseif appendindex == 1
    output.index = index
    out = output;
end


warning('ON','MATLAB:divideByZero');