function out = calclinfielddistance(index, excludeperiods, spikes, linpos, binsize)
%[distance, peak1, peak2] = calclinfielddistance(index, excludeperiods, spikes, linpos, binsize)
%
%Calculates the linear occupancy normalized firing rate for the two cells and
%the distnace between the peaks. The information is not seperated by
%trajectories.
%
%spikes - the 'spikes' cell array for the day you are analyzing
%linpos - the output of LINEARDAYPROCESS for the day you are analyzing. 
%index - [day epoch tetrode cell], two rows (for two cells)
%binsize- the length of each spatial bin (default 2cm)
%excludeperiods - [start end] times for each exlcude period
%


if (nargin < 5)
    binsize = 2;
end

warning('OFF','MATLAB:divideByZero');
statematrix = linpos{index(1,1)}{index(1,2)}.statematrix;
intersections = linpos{index(1,1)}{index(1,2)}.wellSegmentInfo.distanceToIntersection*1000;
timestep = statematrix.time(2,1) - statematrix.time(1,1);
spikes1 = spikes{index(1,1)}{index(1,2)}{index(1,3)}{index(1,4)}.data;
spikes2 = spikes{index(2,1)}{index(2,2)}{index(2,3)}{index(2,4)}.data;

%keep the time and the posindex for each valid spike
spikes1 = spikes1(~isExcluded(spikes1(:,1),excludeperiods),1);
spikes2 = spikes2(~isExcluded(spikes2(:,1),excludeperiods),1);

%bin the distance from two of the wells, where 0 is the W intersection
goodtimes = find(~isExcluded(statematrix.time, excludeperiods));
behavetimes = statematrix.time(goodtimes);
welldist = statematrix.linearDistanceToWells(goodtimes,1:2);   
bins = [];
bins{1} = [fliplr([intersections(1):-binsize:0]) [intersections(1)+binsize:binsize:max(welldist(:,1))]]';
bins{2} = [fliplr([intersections(2):-binsize:0]) [intersections(2)+binsize:binsize:max(welldist(:,2))]]';
welldistbins = [];
for i = 1:2
    welldistbins = [welldistbins bins{i}(lookup(welldist(:,i),bins{i}))];
end

%get all the unique bins
[uniquebins, crap, uniquebinindex] = unique(welldistbins,'rows');

%calculate the distance between each pair of bins
distancematrix = zeros(size(uniquebins,1),size(uniquebins,1));
for i = 1:size(uniquebins,1)
    tmp = repmat(uniquebins(i,:),size(uniquebins,1),1);
    distancematrix(:,i) = max(abs(tmp-uniquebins),[],2);
end

%create a smoothing matrix based on a gaussian curve
distancematrix2 = round(distancematrix/binsize)+1;
smoothcurve = gaussian(2,max(distancematrix2(:))*2);
smoothcurve = smoothcurve(round(length(smoothcurve)/2):end);
smoothcurve(end+1) = 0; 
distancematrix2 = smoothcurve(distancematrix2);
distancematrix2 = distancematrix2 * diag(1./sum(distancematrix2)); %normalize each bins smoothing to 1

%calculate the smoothed occupancy
occupancy = rowcount(uniquebins,welldistbins)*timestep;
smoothoccupancy = distancematrix2'*occupancy;

%calculate the smoothed spike counts
spikecount1 = rowcount(uniquebins,uniquebins(uniquebinindex(lookup(spikes1,behavetimes)),:));
smoothspikecount1 = distancematrix2'*spikecount1;
spikecount2 = rowcount(uniquebins,uniquebins(uniquebinindex(lookup(spikes2,behavetimes)),:));
smoothspikecount2 = distancematrix2'*spikecount2;

%calculate occupancy normalized firing rates
occnormrate1 = smoothspikecount1./smoothoccupancy;
occnormrate2 = smoothspikecount2./smoothoccupancy;

%calculate the distance bateen the peaks
[peak1, peakind1] = max(occnormrate1);
[peak2, peakind2] = max(occnormrate1);
peakdistance = distancematrix(peakind2,peakind1);

out = [peakdistance peak1 peak2];

