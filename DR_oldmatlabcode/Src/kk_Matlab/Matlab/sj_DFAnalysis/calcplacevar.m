function out = calcplacevar(ind, excludetimes, spikes, linpos, cellinfo, bin, minspikes)
%function out = calcplacevar(indices, excludetimes, spikes, linpos, cellinfo, %bin, minspikes)
% Calculates the residuals of the number of spikes in time bins of size bin (in
% seconds) where the residual is the difference between the predicted number 
% of spikes from the linear place field and the actual number fired.
% Only cells with >= minspikes non-excluded spikes are included
% This function then calculates the correlations for all cell pairs in the
% index list, sorted by cell location
%

% get a list of times for the integration of the rate function
times = linpos{ind(1,1)}{ind(1,2)}.statematrix.time;
ntimes = length(times);

% get the list of trajectories the animal was on
traj = linpos{ind(1,1)}{ind(1,2)}.statematrix.traj;

% we need to add the times associated with invalid trajectories to the exclude
% time list.
addexclude = traj >= 1;
tmpe = getExcludePeriods(times, addexclude);
excludetimes = combineExcludePeriods(tmpe, excludetimes);

% get the linear distance at each of those points
lindist = linpos{ind(1,1)}{ind(1,2)}.statematrix.lindist;

% get the time indeces corresponding to the bin edges. 
bins = times(1):bin:times(end);
tbinind = lookup(bins, times);

ncells = size(ind,1);


tmpi = 1;
out.cellloc = {};
for i = 1:ncells
    % check to see if there are enough non-excluded spikes 
    if (length(find(~isExcluded(spikes{ind(i,1)}{ind(i,2)}{ind(i,3)}{ind(i,4)}.data(:,1), excludetimes))) >= minspikes);
	% get the cell locations so we can sort the index list
	out.cellloc{tmpi} = cellinfo{ind(i,1)}{ind(i,2)}{ind(i,3)}{ind(i,4)}.area;
	newind(tmpi,:) = ind(i,:);
	tmpi = tmpi+1;
    end
end
if (~isempty(out.cellloc))
    [out.cellloc cind] = sort(out.cellloc);

    %rearrange the indeces to put the cells in order by location
    ind = newind(cind,:);
    ncells = length(cind);
else
    ncells = 0;
end

resid = zeros(ncells, length(bins) - 1);
rate = zeros(ntimes, 1);
for i = 1:ncells
    % Calculate the linear place fields
    lf{i} = filtercalclinfields(ind(i,1:4), excludetimes, spikes, linpos);

    ntraj = length(lf{i}.trajdata);
    % get rid of any trajectories for which we don't have the linear fields
    for t = 1:max(traj)
        if ((t > ntraj) | isempty(lf{i}.trajdata{t}))
            traj(find(traj == t)) = 0;
        end
    end


    lindistind = zeros(ntraj, ntimes);

    % get, for each trajectory, the rate for each distance
    for t = 1:ntraj
        try
        	lindistrate(t,:) = lf{i}.trajdata{t}(lookup(lindist, lf{i}.trajdata{t}(:,1)),5);
        end
    end

    % get the index, at each time, into the lindistrate matrix
    validtraj = find(traj > 0);
    validtimeind = 1:ntimes;
    validtimeind = validtimeind(validtraj)';
    tmpind = sub2ind(size(lindistrate), traj(validtraj), validtimeind);

    % create a firing rate function across time
    rate(validtraj) = lindistrate(tmpind);

    % add exclude times for all NaN rates
    tmpe = getExcludePeriods(times, isfinite(rate));
    tmpexclude = combineExcludePeriods(excludetimes, tmpe); 

    % get rid of rates associated with excluded times (set to 0 so they don't
    % contribute to the integral
    excluderate = isExcluded(times, tmpexclude);
    rate(find(excluderate)) = 0;

    % get the cumulative integral the rate function 
    rateint = cumtrapz(times, rate);

    % the number of spikes is the difference in the cumulative integral from
    % the beginning to the end of each bin.
    nexpected = diff(rateint(tbinind));

    % get the included spikes
    s = spikes{ind(i,1)}{ind(i,2)}{ind(i,3)}{ind(i,4)};
    s.data = s.data(~isExcluded(s.data(:,1), tmpexclude), 1);

    % histogram the spikes to get the actual number in each bin
    [nspikes tmp] = flhist(s.data, bins);
    resid(i,:) = nexpected - nspikes';
end
out.resid = resid;

% create an matrix with all pairwise correlations
out.corr = zeros(ncells);
out.corrp = zeros(ncells);
if (exist('lf'))
    out.lf = lf
else
    out.lf = {};
end
for i = 1:ncells
    for j = i+1:ncells
	% find the elements of the residuals that are both nonzero
	nonzero = find(resid(i,:) & resid(j,:));
	if (length(nonzero) > 5)
	    [out.corr(i,j) out.corrp(i,j)] = corr(resid(i,nonzero)', ...
					resid(j,nonzero)', 'type', 'Spearman');
            %if (out.corrp(i,j) < 0.01)
	%	plot(resid(i,nonzero));
	%	hold on
	%	plot(resid(j,nonzero));

	else
	    out.corr(i,j) = NaN;
	    out.corrp(i,j) = NaN;
	end
    end
end





