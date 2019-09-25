function out = calcresiduals_evensampling(ind, excludetimes, spikes, linpos, pos, bin, minspikes)
%function out = calcresiduals_evensampling(indices, excludetimes, spikes, linpos, pos, bin, minspikes)
% Calculates the residuals of the number of spikes in time bins of size bin
% (in seconds) where the residual is the difference between the predicted 
% number of spikes from the linear place field and the actual number fired.
% Only cells with >= minspikes for non-excluded times are included
% This function then calculates the correlations for all cell pairs in the
% index list, sorted by cell location
% This function also computes the mean speed for each time bin.

% This function is the same as calcresiduals EXCEPT that it restricts the
% calculation to locations where speed was sampled at multiple speeds

% get a list of times for the integration of the rate function
times = linpos{ind(1,1)}{ind(1,2)}.statematrix.time;
ntimes = length(times);

% get the list of trajectories the animal was on
traj = linpos{ind(1,1)}{ind(1,2)}.statematrix.traj{6};

% add the times associated with invalid trajectories to the exclude times
addexclude = traj >= 1;
tmpe = getExcludePeriods(times, addexclude);
excludetimes = combineExcludePeriods(tmpe, excludetimes);

% get the linear distance
lindist = linpos{ind(1,1)}{ind(1,2)}.statematrix.lindist;
distbin = 0:2.5:max(lindist);

% get the time indeces corresponding to the bin edges. 
bins = times(1):bin:times(end);
tbinind = lookup(bins, times);

% get the speed for each time bin
vel = pos{ind(1,1)}{ind(1,2)}.data(lookup(times,pos{ind(1,1)}{ind(1,2)}.data(:,1)),8);
traj_dist = zeros(size(traj));
traj_dist(traj == 1| traj==2) = 1;
traj_dist(traj == 3| traj==4) = 2;
dist = nan(2,length(tbinind)-1);
speed = zeros(length(tbinind)-1,1);
traj_distance = zeros(size(speed));
for i = 1:(length(tbinind)-1)
    speed(i) = mean(vel(tbinind(i):tbinind(i+1)));
    traj_distance(i) = round(mean(traj_dist(tbinind(i):tbinind(i+1))));
    if traj_distance(i)~=0
        dist(traj_distance(i),i) = lookup(mean(lindist(tbinind(i):tbinind(i+1))),distbin);
    end
end

% Figure out which locations were sampled at multiple speeds
speedbin = ones(size(speed));
speedbin(speed > 4) = 2;
validdist = hist(dist(:,speedbin==1)',1:length(distbin))>0 & hist(dist(:,speedbin==2)',1:length(distbin))>0;

valid = zeros(size(speedbin));
for i = 1:length(dist)
    %Determine if the distance is defined at this time
    if traj_distance(i) == 1 || traj_distance(i) == 2
        %Determin if this distance is valid
        valid(i) = validdist(dist(traj_distance(i),i),traj_distance(i));
    end
end
valid = ones(size(speedbin));

ncells = size(ind,1);   
% get the cells that have the minimum number of spikes
valid_cells = zeros(size(ncells));
for i = 1:ncells
    if ~isempty(spikes{ind(i,1)}{ind(i,2)}{ind(i,3)}{ind(i,4)}.data)
        if (length(find(~isExcluded(spikes{ind(i,1)}{ind(i,2)}{ind(i,3)}{ind(i,4)}.data(:,1), excludetimes))) >= minspikes);
            valid_cells(i) = 1;
        end
    end
end
ind = ind(logical(valid_cells),:);
ncells = sum(valid_cells);

%initialize variables
resid = zeros(ncells, length(bins) - 1);
peak_rate = zeros(size(ncells));
lf = cell(size(ncells));

for i = 1:ncells
    rate = zeros(ntimes, 1);
    % Calculate the linear place fields
    lf{i} = filtercalclinfields(ind(i,1:4), excludetimes, spikes, linpos);

    ntraj = length(lf{i}.trajdata);
    % get rid of any trajectories for which we don't have the linear fields
    for t = 1:max(traj)
        if ((t > ntraj) || isempty(lf{i}.trajdata{t}))
            traj(traj == t) = 0;
        end
    end

    lindistrate = zeros(ntraj, ntimes);
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
    peak_rate(i) = max(rate);
    
    % add exclude times for all NaN rates
    tmpe = getExcludePeriods(times, isfinite(rate));
    tmpexclude = combineExcludePeriods(excludetimes, tmpe); 

    % get rid of rates associated with excluded times (set to 0 so they don't
    % contribute to the integral
    excluderate = isExcluded(times, tmpexclude);
    rate(logical(excluderate)) = 0;

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

if ncells >= 2
    pairsind = nchoosek(1:ncells,2);
    tmpout = nan(size(pairsind,1),9);
    for i = 1:size(pairsind,1)
        % Column 1 & 2: Cell1 Cell2 (index into cindex)
        tmpout(i,[1 2]) = pairsind(i,:);

        % Column3: Cell1 peak rate
        tmpout(i,3) = peak_rate(pairsind(i,1));

        % Column4: Cell2 peak rate
        tmpout(i,4) = peak_rate(pairsind(i,2));

        % Column5: Cell pair normalized overlap
        tmpout(i,5) = calcoverlap(lf{pairsind(i,1)}.trajdata,lf{pairsind(i,2)}.trajdata);

        resid1 = resid(pairsind(i,1),:);
        resid2 = resid(pairsind(i,2),:);
        nonzero = (resid1 & resid2)';
        if sum(~isnan(resid1(nonzero)) & ~isnan(resid2(nonzero)))*bin >10
            % Column 6 & 7: correlation coefficient and pvalue
            [tmpout(i,6) tmpout(i,7)] = corr(resid1(nonzero)',resid2(nonzero)', ...
                'type','Spearman','rows','complete');
            for s = 1:max(speedbin)
                if sum(~isnan(resid1(valid & nonzero & speedbin == s)) &...
                        ~isnan(resid2(valid & nonzero & speedbin == s)))*bin > 10
                    % Column 8 - 9: speed correlation coeffiecients
                    tmpout(i,7+s) = corr(resid1(valid & nonzero & speedbin == s)',...
                        resid2(valid & nonzero & speedbin == s)',...
                        'type','Spearman','rows','complete');
                end
            end
        end
    end                        
elseif ncells == 1
    tmpout = nan(1,9);
    tmpout([1 3]) = [1 peak_rate(1)];
else
    tmpout = [];
end 
    
out.speed = speed;
out.cindex = ind;
out.resid = resid;
out.lf = lf;
out.cell_pair_info = tmpout;