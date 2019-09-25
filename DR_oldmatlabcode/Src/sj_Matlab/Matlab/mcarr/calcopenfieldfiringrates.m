function output = calcopenfieldfiringrates(index, excludetimes, spikes,pos, task, varargin)
%[out] = openfieldoccupancy(index, excludetimes, spikes,pos, options)
%
%Calculates the 2d occupancy normalized firing rate for the cell.  Bases
%the bins on the linear coordinates provided in task structure.
%
% options:
%   binsize- the length of each spatial bin (default 1 cm)
%   std - defines the shape of the 2d gaussian used to smooth spikerate.
%              (default 1)
%   appendindex - 0 or 1, 1 appends index infront of output (default 1 )
%
%The output is a structure with fields: occupancy, bin vector x (.xticks), bin vector
%y (.yticks), bin spike count (.spikes), occ normailized firing per bin (.spikerate), and smoothed occ
% normalized firing (.smoothedspikerate).


appendindex = 1;

output.type= nan(size(index,1),1);
output.quadrants = nan(size(index,1),1);
output.peak = nan(size(index,1),1);
output.rate = cell(size(index,1),1);
output.meanrate = nan(size(index,1),1);
output.quadmeanrate = nan(size(index,1),4);
output.location = cell(size(index,1),1);
output.quadmapping = nan(1,4);
binsize = 5;

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

%Define variables and filter out exclude times
pos = pos{index(1,1)}{index(1,2)}.data;
task = task{index(1,1)}{index(1,2)};
valid =  ~isExcluded(pos(:,1), excludetimes);

timestep = median(diff(pos(:,1)));
totalontime = sum(valid)*timestep;
quadtypes = task.quadrants;
pos = pos(valid,:);

type = zeros(4,1);
%Determine the type of each quadrant
for q = 1:4
    %Determine what type of quadrant q is
    if task.objects(q) && task.novel(q);
        type(q) = 1;
        quadtypes(task.quadrants==q) = 1;
    elseif task.objects(q) && ~task.novel(q);
        type(q) = 2;
        quadtypes(task.quadrants==q) = 2;
    elseif ~task.objects(q) && task.novel(q);
        type(q) = 3;
        quadtypes(task.quadrants==q) = 3;
    elseif ~task.objects(q) && ~task.novel(q);
        type(q) = 4;
        quadtypes(task.quadrants==q) = 4;
    else
        type(q) = NaN;
    end   
end
quadtime = hist(quadtypes(valid),[1 2 3 4])*timestep;

%Calcualte preference, if this is a familiar session, compare two objects
if any(1 == type)
    obj = [find(type==1) find(type==2)];
    output.preference = (sum(task.quadrants(valid)==obj(1))-sum(task.quadrants(valid)==obj(2)))./(sum(task.quadrants(valid)==obj(1))+sum(task.quadrants(valid)==obj(2)));
else
    obj = find(type ==2);
    output.preference = (sum(task.quadrants(valid)==obj(1))-sum(task.quadrants(valid)==obj(2)))./(sum(task.quadrants(valid)==obj(1))+sum(task.quadrants(valid)==obj(2)));
end

%Define spatial bins
bincenters = task.linearcoords;
minxy = min(pos(:,2:3));
maxxy = max(pos(:,2:3));
binx = bincenters(1)-ceil((bincenters(1)-minxy(1))/binsize)*binsize:binsize:bincenters(1)+ceil((-bincenters(1)+maxxy(1))/binsize)*binsize;
biny = bincenters(2)-ceil((bincenters(2)-minxy(2))/binsize)*binsize:binsize:bincenters(2)+ceil((-bincenters(2)+maxxy(2))/binsize)*binsize;
edges{1} = binx; edges{2} = biny;
    
occupancy = hist3(pos(:,2:3),'Edges',edges);
occupancy = occupancy.*timestep;
invalid = occupancy < 0.25;
occupancy(invalid) = NaN;
output.xticks = binx;
output.yticks = biny;

%Go through each cell and compute the firing rate
for cells = 1:size(index,1)
        
    s = spikes{index(cells,1)}{index(cells,2)}{index(cells,3)}{index(cells,4)}.data;
    
	indgoodspikes = ~isExcluded(s(:,1), excludetimes);
    spiketimes = s(indgoodspikes,1);
    spikequads = s(indgoodspikes,7);
    if ~isempty(spiketimes)
        
        %Histogram the spikes and compute firing rate
        tmpspikes =hist3(s(indgoodspikes,2:3),'Edges',edges);
        output.rate{cells} = tmpspikes./occupancy;

        output.peak(cells) = max(max(output.rate{cells}));
        
        %Figure out the center of mass of the place field and quadrant type
        val = output.rate{cells}(~isnan(output.rate{cells}));
        
        [x y] = find(~isnan(output.rate{cells}));
        sx = sum(val.*binx(x)')./sum(val);
        sy = sum(val.*biny(y)')./sum(val);
        %Determine the COM location relative to each quadrant center
        quadcenters = [median(binx(binx<bincenters(1))) median(binx(binx>bincenters(1))) ...
            median(biny(biny<bincenters(1))) median(biny(biny>bincenters(2)))];
        output.location{cells}(1,:) = [sx-quadcenters(1) sy-quadcenters(4)];
        output.location{cells}(2,:) = [sx-quadcenters(2) sy-quadcenters(4)];
        output.location{cells}(3,:) = [sx-quadcenters(2) sy-quadcenters(3)];
        output.location{cells}(4,:) = [sx-quadcenters(1) sy-quadcenters(3)];
        
        if any(isnan(output.location{cells}))
            keyboard
        end
        
        if sx <= bincenters(1) && sy >= bincenters(2)
            output.type(cells) = type(1);
            output.quadrants(cells) = 1;
        elseif sx>bincenters(1) && sy >= bincenters(2)
            output.type(cells) = type(2);
            output.quadrants(cells) = 2;
        elseif sx>bincenters(1) && sy < bincenters(2)
            output.type(cells) = type(3);
            output.quadrants(cells) = 3;
        elseif sx<= bincenters(1) && sy < bincenters(2)
            output.type(cells) = type(4);
            output.quadrants(cells) = 4;
        end
        
        output.quadmapping = type;
        %Determine mean firing rate
        output.meanrate(cells) = length(spiketimes)./totalontime;
        
        %Determine mean firing rate by quadrant
        spikequads = hist(quadtypes(spikequads),[1 2 3 4]);
        output.quadmeanrate(cells,:) = spikequads./quadtime;
    end

end

%Reshape the output.quadmeanrate so that the order is [1 2 3 4] or [2 2 4 4]

if appendindex
    output.index = index;
end

if isfield(task,'familiarsession')
    output.familiarsession = task.familiarsession;
else
    output.familiarsession = [];
end

if isfield(task,'novelsession')
    output.novelsession = task.novelsession;
else
    output.novelsession = [];
end

end