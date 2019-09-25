function output = calcmeanrates_by_quadrant(index, excludetimes, spikes,pos, task, varargin)
%[out] = openfieldoccupancy(index, excludetimes, spikes,pos, options)
%
%Calculates the population mean firing rate by quadrant
%
% options:
%   appendindex - 0 or 1, 1 appends index infront of output (default 1 )
%
%The output is a structure with fields: occupancy, bin vector x (.xticks), bin vector
%y (.yticks), bin spike count (.spikes), occ normailized firing per bin (.spikerate), and smoothed occ
% normalized firing (.smoothedspikerate).


appendindex = 1;

output.type= nan(size(index,1),1);
output.quadmeanrate = nan(size(index,1),4);
output.totalmeanrate = nan(size(index,1),1);

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

%Define variables and filter out exclude times
pos = pos{index(1,1)}{index(1,2)}.data;
task = task{index(1,1)}{index(1,2)};
valid =  ~isExcluded(pos(:,1), excludetimes);

timestep = median(diff(pos(:,1)));
quadtime = hist(task.quadrants(valid),[1 2 3 4]).*timestep;
totaltime = sum(valid).*timestep;
pos = pos(valid,:);
type = zeros(4,1);
%Determine the type of each quadrant
for q = 1:4
    %Determine what type of quadrant q is
    if task.objects(q) && task.novel(q);
        type(q) = 1;
    elseif task.objects(q) && ~task.novel(q);
        type(q) = 2;
    elseif ~task.objects(q) && task.novel(q);
        type(q) = 3;
    elseif ~task.objects(q) && ~task.novel(q);
        type(q) = 4;
    else
        type(q) = NaN;
    end   
end

tmpout = zeros(size(index,1),5);
%Go through each cell and gather the spikes by quadrant
for cells = 1:size(index,1)
        
    s = spikes{index(cells,1)}{index(cells,2)}{index(cells,3)}{index(cells,4)}.data;
    
	indgoodspikes = ~isExcluded(s(:,1), excludetimes);
    spiketimes = s(indgoodspikes,7);
    if ~isempty(spiketimes)
        
        %Histogram the spikes and compute firing rate
        tmpspikes = hist(task.quadrants(spiketimes,1),[1 2 3 4]);

        tmpout(cells,:) = [tmpspikes./quadtime length(spiketimes)./totaltime];
    end

end

output.quadmeanrate = tmpout(:,1:4);
output.totalmeanrate = tmpout(:,5);
output.type = type;

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