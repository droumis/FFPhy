function out = calcripplespikingvelocity(index, excludeperiods, spikes, ripples, pos, cellinfo, varargin)

%Set options
cellfilter = [];

%Process options
for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})
        switch(varargin{option})
            case 'cellfilter'
                cellfilter = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

%Find valid ripples
if isempty(cellfilter)
    riptimes = getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
        '(isequal($area, ''CA1''))',...
        'excludeperiods', excludeperiods,'minstd',3);
    tet = evaluatefilter(cellinfo{index(1,1)}{index(1,2)},'isequal($area,''CA1'')');
    tet = unique(tet(:,1));
else
    riptimes = getripples([index(1,1) index(1,2)], ripples, cellinfo, 'cellfilter', ...
        cellfilter,'excludeperiods', excludeperiods,'minstd',3);
    tet = evaluatefilter(cellinfo,'cellfilter');
    tet = unique(tet(:,3));
end

%Remove ripples that are too close to the beginning or end
while riptimes(1,3)<pos{index(1,1)}{index(1,2)}.data(1,1) + 0.5
	riptimes(1,:) = [];
end
while riptimes(end,3)> pos{index(1,1)}{index(1,2)}.data(end,1)- 0.5
    riptimes(end,:) = [];
end

%Find velocity for each valid ripple
speed = nan(size(riptimes(:,1)));
time_start = lookup(riptimes(:,1),pos{index(1,1)}{index(1,2)}.data(:,1));
time_end = lookup(riptimes(:,2),pos{index(1,1)}{index(1,2)}.data(:,1));
for s = 1:length(time_start)
    speed(s) = mean(pos{index(1,1)}{index(1,2)}.data(time_start(s):time_end(s),8));
end

bin = [1/2 1 2 4 8 16 32];
speedbin = lookup(speed,bin);
invalid = hist(speedbin,1:length(bin))<5;
spike_matrix = zeros(size(index,1),size(riptimes,1));
prop_active = nan(size(index,1),length(bin));
total_prop_active = nan(size(index,1),1);

%Go through each cell and determine its firing rate during each ripple
for cellcount = 1:size(index,1)
    ind = index(cellcount,:);
    if ~isempty(spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data)
        spiketimes = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data(:,1);
    else
        spiketimes = [];
    end
    %Find modulation during ripples
    valid = find(~isExcluded(spiketimes, excludeperiods));
    spikebins = periodAssign(spiketimes(valid), riptimes(:,[1 2]));
    rate = hist(spikebins,0:1:size(riptimes,1));
    spike_matrix(cellcount,:) = rate(2:end)'./(riptimes(:,2)-riptimes(:,1));
    prop_active(cellcount,:) = accumarray(speedbin,rate(2:end)>0,[length(bin) 1], @(x) sum(x)./length(x));
    total_prop_active(cellcount) = sum(rate(2:end)>0)./size(riptimes,1);
end
prop_active(:,invalid) = NaN;

%Go through each valid tetrode and determine its maxthresh

maxthresh_matrix = nan(length(tet),size(riptimes,1));
maxthresh_total = nan(length(tet),size(riptimes,1));

for tetcount = 1:length(tet)
    ripple_ind = lookup(ripples{index(1,1)}{index(1,2)}{tet(tetcount)}.starttime,riptimes(:,1));
    maxthresh_matrix(tetcount,ripple_ind) = ripples{index(1,1)}{index(1,2)}{tet(tetcount)}.maxthresh;
    std_epoch = ripples{index(1,1)}{index(1,2)}{tet(tetcount)}.std;
    baseline_epoch = ripples{index(1,1)}{index(1,2)}{tet(tetcount)}.baseline;
    tmp = maxthresh_matrix(tetcount,:).*std_epoch + baseline_epoch;
    
    std_day = ripples{index(1,1)}{index(1,2)}{tet(tetcount)}.totalstd;
    baseline_day = ripples{index(1,1)}{index(1,2)}{tet(tetcount)}.totalbaseline;
    maxthresh_total(tetcount,:) = (tmp - baseline_day)./std_day;
end


out.maxthresh = nanmean(maxthresh_matrix,1);
out.totalmaxthresh = nanmean(maxthresh_total,1);
out.speed = speed;
out.active_prob = prop_active;
out.total_proportion_active = total_prop_active;
end

function out = periodAssign(times, periods)
% out = periodAssign(times, periods)
% TIMES is a vector of times
% PERIODS is an N by 2 list of start and end times
% Returns the index of the period that each time falls into.  If a time
% does not fall into one of the periods, a zero is returned.
% This function assumes that the periods are not overlapping.
%

if ~isempty(periods)
    oneborder = [(periods(:,1)-.0000001);periods(:,2)+.0000001];
    oneborder(:,2) = 0;
    insideborder = [(periods(:,1)+.0000001) (1:length(periods))'; ...
        (periods(:,2)-.0000001) (1:length(periods))'];    
    sortedMatrix = [[-inf 0]; sortrows([oneborder;insideborder],1); ...
        [inf 0]];
else
    sortedMatrix = [[-inf 0]; [inf 0]];
end
out = sortedMatrix(lookup(times,sortedMatrix(:,1)),2);
end