function out = calcspatialrate(indices, excludetimes, spikes, pos, varargin)


appendindex = 0;
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end




for cellNum = 1:size(indices.spikes,1)
    index = [indices.epochs indices.spikes(cellNum,:)];
    out(cellNum).xloc = [];
    out(cellNum).yloc = [];
    out(cellNum).occupancy = []; 
    out(cellNum).spikes = [];
    out(cellNum).spikerate = [];
    out(cellNum).smoothedspikerate = [];    
    out(cellNum).index = [];
    
    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.time)
            spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.time(:,1);
            goodspikes = ((~isExcluded(spiketimes, excludetimes)));
            goodspikedata = spikes{index(1)}{index(2)}{index(3)}{index(4)}.time(find(goodspikes),:);              
            postimes = pos{index(1)}{index(2)}.time(:,1);       
            spikepos = [pos{index(1)}{index(2)}.x(lookup(goodspikedata,postimes)) pos{index(1)}{index(2)}.y(lookup(goodspikedata,postimes))];
            goodspikedata = [goodspikedata spikepos];
            
            goodpostimes = ((~isExcluded(postimes, excludetimes)));
            goodposdata = [pos{index(1)}{index(2)}.time(find(goodpostimes)) pos{index(1)}{index(2)}.x(find(goodpostimes)) pos{index(1)}{index(2)}.y(find(goodpostimes))];
            
            
            
            tmpresult = calc2drateINT(goodspikedata, goodposdata);
            if ~isempty(tmpresult)
                tmpresult.index = index;
                out(cellNum) = tmpresult;
            end
        else
            goodspikes = 0;
            goodspikedata = [];
        end
    else
        goodspikes = nan;
        goodspikedata = [];
    end
    
    
    
end
%-------------------------------------------------------------

function [output] = calc2drateINT(spikes, pos, varargin)
%[output] = twodoccupancy(spikes, pos, index, binsize)
%[output] = twodoccupancy(spikes, pos, index)
%
%Calculates the 2d occupancy normalized firing rate for the cell and
%organizes the output into the different trajectories.
%
%spikes - the data field from the spikes cell array
%pos - the data field from the pos cell array
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

std = 1;
filter = [];
binsize = 1;
goodspikeind = ones(size(spikes,1),1);
goodposind = ones(size(pos,1),1);
output = [];

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'binsize'
                binsize = varargin{option+1};
            case 'std'
                std = varargin{option+1};
            case 'filter'
                filter = varargin{option+1};
            case 'goodspikeind'
                goodspikeind = varargin{option+1};
            case 'goodposind'
                goodposind = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end




timestep = min(diff(pos(:,1))); 

allpos = pos(:,[2 3]);
tmpposition = pos(find(goodposind),[2 3]);
tmpspikes = spikes(find(goodspikeind),[2 3]);

if ~isempty(tmpposition)
    minx = floor(min(allpos(:,1)));
    maxx = ceil(max(allpos(:,1)));
    binx = (minx:binsize:maxx);
    miny = floor(min(allpos(:,2)));
    maxy = ceil(max(allpos(:,2)));
    biny = (miny:binsize:maxy);
    
    
    
   
    
    
    if ~isempty(tmpspikes)
        output.xloc = binx;
        output.yloc = biny;
        
        output.occupancy = hist2d(tmpposition, binx, biny);
        output.spikes = [];
        output.spikerate = [];
        output.smoothedspikerate = [];
        
        [output.spikes] = hist2d(tmpspikes, binx, biny);
        nonzero = find(output.occupancy ~= 0);
        output.spikerate = zeros(size(output.spikes));
        output.spikerate(nonzero) = output.spikes(nonzero) ./(timestep* output.occupancy(nonzero) );
        
        %smoothed occupancy
        g = gaussian2(std,round(6*std));
        output.smoothedspikerate = filter2(g,(output.spikerate));
        smoothedoccupancy = [];
        smoothedoccupancy = zeros(size(output.spikes));
        smoothedoccupancy = filter2(g, output.occupancy);
        zero = find(smoothedoccupancy == 0);
        output.smoothedspikerate(zero) = -1;
        %rate = output.smoothedspikerate;
    else
        %rate = [];
    end
else
    %rate = [];
end
warning('ON','MATLAB:divideByZero');

