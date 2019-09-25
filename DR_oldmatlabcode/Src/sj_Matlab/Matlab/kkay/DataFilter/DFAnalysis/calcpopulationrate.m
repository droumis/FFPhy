function out = calcpopulationrate(indices, excludetimes, spikes, timebin)

if (nargin < 4)
    timebin = 60;
end

epochind = indices(:,1:2);
epochind = unique(epochind,'rows');
if (size(epochind,1) ~= 1)
    error('Indices must have only one unique day/epoch');
end

timerange = spikes{indices(1,1)}{indices(1,2)}{indices(1,3)}{indices(1,4)}.timerange/10000;
timeperiods = [timerange(1):timebin:timerange(2)];
%for each time bin, calculate the total non-excluded time
totaltimes = calcSegmentTimeOn(timeperiods, excludetimes);

numcells = 0;
out = [];
for i = 1:size(indices,1)
    if ~isempty(spikes{indices(i,1)}{indices(i,2)}{indices(i,3)}{indices(i,4)})
        numcells = numcells+1;
        if ~isempty(spikes{indices(i,1)}{indices(i,2)}{indices(i,3)}{indices(i,4)}.data)
            tmpspiketimes = spikes{indices(i,1)}{indices(i,2)}{indices(i,3)}{indices(i,4)}.data(:,1);
            tmpspiketimes = tmpspiketimes(find(~isExcluded(tmpspiketimes, excludetimes)));
        else
            tmpspiketimes = [];
        end
        spikecount = [];
        for i = 1:length(timeperiods)-1
            spikecount(i) = sum((tmpspiketimes >= timeperiods(i)) & (tmpspiketimes < timeperiods(i+1)));
        end
        out = [out; (spikecount./totaltimes)];
    end
end
