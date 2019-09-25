function [out] = getcoactiveripples_and_overlap(index, excludetimes, ripples, pos, task, spikes, cellinfo, cellfilter, varargin)
%[out] = getcoactiveripples_and_overlap(index, excludeperiods,ripples, pos, task, spikes, cellinfo, cellfilter)
%Computes the coactivation z-score for each pair of cells and the overlap
%of the place fields

%   index [day epoch]
%   tetlist is a list of tetrodes to include in the analysis
%
%   options:
%   'appendindex' , 1 or 0, default 0
%           set to 1 to append the cell index to the output [day epoch
%           value]
%   'normalize'
%   'difftetrode'
%   'binsize', size of spatial bins used to calculate the firing rate
%   out

% assign the options
appendindex = 0;

difftetrode = 0;
firingrate = cell(size(index(1,1)));
quad = cell(size(index(1,1)));
binsize = 5;

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'appendindex'
            appendindex = varargin{option+1};
        case 'difftetrode'
            difftetrode = varargin{option+1};
        case 'binsize'
            binsize = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

if size(index,1)<2
    out = [];
elseif  ~isempty(evaluatefilter(cellinfo{index(1,1)}{index(1,2)},'isequal($area,''CA1'')'))
    
    pos = pos{index(1,1)}{index(1,2)}.data;
    task = task{index(1,1)}{index(1,2)};
    valid =  ~isExcluded(pos(:,1), excludetimes);

    timestep = median(diff(pos(:,1)));
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

    %Define spatial bins
    bincenters = task.linearcoords;
    minxy = min(pos(:,2:3));
    maxxy = max(pos(:,2:3));
    binx = bincenters(1)-ceil((bincenters(1)-minxy(1))/binsize)*binsize:binsize:bincenters(1)+ceil((-bincenters(1)+maxxy(1))/binsize)*binsize;
    biny = bincenters(2)-ceil((bincenters(2)-minxy(2))/binsize)*binsize:binsize:bincenters(2)+ceil((-bincenters(2)+maxxy(2))/binsize)*binsize;
    edges{1} = binx; edges{2} = biny;

    occupancy = hist3(pos(:,2:3),'Edges',edges);
    occupancy = occupancy.*timestep;
    invalid = occupancy < 0.5;
    occupancy(invalid) = NaN;

    %Get ripples
    r = getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', 'isequal($area, ''CA1'')','excludeperiods', excludetimes,'minstd',3);

    %for each cell compute if active in each ripple
    ripspikes = zeros(size(index,1), size(r,1));
    for c = 1:size(index,1) %for each cells
        if ~isempty(spikes{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}) && ~isempty(spikes{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}.data)
            %get activation in ripples
            spiketimes = spikes{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}.data(:,1);
            indgoodspikes = find(~isExcluded(spiketimes, excludetimes));
            spikebins = periodAssign(spiketimes(indgoodspikes), r(:,[1 2]));

            ripactive = unique(spikebins(spikebins > 0));
            ripspikes(c, ripactive) = 1; %if spike in a bin for this cell, update ripspikes, each row is a cell, each column is a ripple

            %Calculate firing rate
            tmpspikes =hist3(spikes{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}.data(indgoodspikes,2:3),'Edges',edges);
            firingrate{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)} = tmpspikes./occupancy;
            
            %Figure out the center of mass of the place field and quadrant type
            sx = nansum((nanmean(firingrate{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}').*binx)/nansum(nanmean(firingrate{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}')));
            sy = nansum((nanmean(firingrate{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)}).*biny)/nansum(nanmean(firingrate{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)})));
            
            if sx == 0 && sy == 0
                quad{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)} = 0;
            elseif sx <= bincenters(1) && sy >= bincenters(2)
                quad{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)} = type(1);
            elseif sx>bincenters(1) && sy >= bincenters(2)
                quad{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)} = type(2);
            elseif sx>bincenters(1) && sy < bincenters(2)
            	quad{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)} = type(3);
            elseif sx<= bincenters(1) && sy < bincenters(2)
                quad{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)} = type(4);
            end
        else
            firingrate{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)} = [];
            quad{index(c,1)}{index(c,2)}{index(c,3)}{index(c,4)} = NaN;
        end
        
    end

    %for each pair figure out how often coactive in ripples & overlap
    pairs = nchoosek(1:size(index,1), 2);
    if difftetrode == 1
    	difftet = logical(index(pairs(:,1),3) - index(pairs(:,2),3)); %identify cells pairs with cells on different tetrode
        pairs = pairs(difftet,:); %only include cells pairs with cells on different tetrode
    end
    coactivrip = zeros(size(pairs,1),1); overlap = coactivrip; coactz = coactivrip; jtsurprise = coactivrip;
    quadrants = zeros(size(pairs,1),2);
	for q = 1:size(pairs,1)  %for each pair figure out how often coactive in ripples and in run
        actAB = length(find(sum(ripspikes(pairs(q,:),:),1)==2)) ; %number of ripples where both cells active
        coactivrip(q) = actAB/size(ripspikes,2); % divided by number of ripples
        actA = sum(ripspikes(pairs(q,1),:)>=1);
        actB = sum(ripspikes(pairs(q,2),:)>=1);
        N = size(ripspikes,2);
        coactz(q) = coactivezscore(N, actAB, actA, actB);
        jtsurprise(q) = jointsurprise(N, actAB, actA, actB);
        
        %Calculate overlap
        if ~isempty(firingrate{index(1,1)}{index(1,2)}{index(pairs(q,1),3)}{index(pairs(q,1),4)}) && ~isempty(firingrate{index(1,1)}{index(1,2)}{index(pairs(q,2),3)}{index(pairs(q,2),4)})
            rate1 = firingrate{index(1,1)}{index(1,2)}{index(pairs(q,1),3)}{index(pairs(q,1),4)};
            rate2 = firingrate{index(1,1)}{index(1,2)}{index(pairs(q,2),3)}{index(pairs(q,2),4)};
        
            overlap(q) = calc2doverlap(rate1,rate2, 'binsize', binsize);
            quadrants(q,:) = [quad{index(1,1)}{index(1,2)}{index(pairs(q,1),3)}{index(pairs(q,1),4)} quad{index(1,1)}{index(1,2)}{index(pairs(q,2),3)}{index(pairs(q,2),4)}];
        end
    end
        
	if appendindex == 0
        out = [overlap coactivrip coactz repmat(size(r,1), length(overlap), 1) jtsurprise quadrants];
	elseif appendindex ==1
        out = [repmat(index(1,1:2), length(overlap), 1) overlap coactivrip coactz repmat(size(r,1), length(overlap), 1) jtsurprise quadrants];
    end
else
    out = [];
end



