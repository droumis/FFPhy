function trajdata = sj_calclinfields(spikes,statevector,lindist,linpos,pos,index,binsize, minabsvel)

% Shantanu: Adding Abs Vel from "pos" field. Accurate calculation of
% velocity, and independent of minvel inputed/used in sj_getbehavestate
% 10 Jun 2011

%trajdata = CALCLINFIELDS(spikes,statevector,lindist,linpos, index)
%trajdata = CALCLINFIELDS(spikes,statevector,lindist,linpos,index,binsize)
%
%Calculates the linear occupancy normalized firing rate for the cell and
%organizes the output into the different trajectories.
%
%spikes - the 'spikes' cell array for the day you are analyzing
%statevector - the outputs of GETBEHAVESTATE. This is a vector with the traj
%              number for each position time (1 based). -1 values signify
%              invalid times and are not used.
%lindist     - also an output of GETBEHAVESTATE. Linear distance traveled.
%linpos - the output of LINEARDAYPROCESS for the day you are analyzing.
% pos data for day in 'pos' file
%index - [day epoch tetrode cell]
%binsize- the length of each spatial bin (default 2cm)
% minabsvel - for runs and place fields

%The output is a cell array where each cell contains a matrix
%descibing one trajectory.  These matrices are n by 5, where n is the
%number of spatial bins. The columns are: 
% 1)linear bin location,
% 2) bin occupancy (seconds),
% 3) bin spike count, 
% 4) occ normailized firing per bin, and
% 5) smoothed occ normalized firing. 
%If the cell is empty, the animal did not enter that trajectoryor no spikes for that trajectory.
%

if (nargin < 7)
    binsize = 2;
    %numspatialbins = 100;
else
    %numspatialbins = binsPerTraj;
end

if nargin < 8,
    minabsvel = 3; % cm/sec - Most conservative for runs and place fields
end

% Use Calebs speed filter if necessary
defaultfilter = 'velocitydayprocess_filter.mat';
eval(['load ', defaultfilter]);

warning('OFF','MATLAB:divideByZero');
statematrix = linpos{index(1)}{index(2)}.statematrix;
wellSegmentInfo = linpos{index(1)}{index(2)}.wellSegmentInfo;
segmentInfo = linpos{index(1)}{index(2)}.segmentInfo;
trajwells = linpos{index(1)}{index(2)}.trajwells;
spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
spikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;

% Position Data with Fields 'time x y dir vel x-sm y-sm dir-sm vel-sm'
posdata = pos{index(1)}{index(2)}.data;
if size(posdata,2)>5 % already smoothed position and filtered velocity
    abspos = posdata(:,6:7);
    absvel = abs(posdata(:,9)); % Should already be absolute value
else
    abspos = posdata(:,2:3);
    currvel = abs(pos{day}{epoch}.data(:,5));
    % If vel-sm does not exist (if forgotten), create it
    absvel = [filtfilt(velocityfilter.kernel, 1, currvel)];
    %smooth_currvel = smoothvect(currvel,filt); 2nd option for convolution
end
% Using original velocity
absvel = abs(posdata(:,5));

if length(statevector)~=length(absvel)
    finalind = min([length(statevector) length(absvel)]);
    statevector = statevector(1:finalind);
    absvel = absvel(1:finalind);
    abspos = abspos(1:finalind);
end

if (length(segmentInfo.segmentLength) > 1)
    for i = 1:size(trajwells,1)
        %calculate the linear length of the trajectory
        trajlength(i) = sum(segmentInfo.segmentLength(wellSegmentInfo.pathTable{trajwells(i,1),wellSegmentInfo.segmentIndex(trajwells(i,2))}));
    end
else
    trajlength(1) = segmentInfo.segmentLength(1);
end


%posindexfield = isdatafield(spikesfields,'posindex');
posindexfield = 7;
if ~isempty(spikes)
    spikes = spikes(:,[1 posindexfield]);
    spikes(:,3) = statevector(spikes(:,2)); %add the traj state for each spike
    spikes(:,4) = lindist(spikes(:,2)); %add the linear distance
    spikes(:,5) = absvel(spikes(:,2)); %add the absolute velocity - Should be the same size as statevector
else
    spikes = [0 0 -1 0];
end

trajnum = max(statevector);
timestep = statematrix.time(2,1) - statematrix.time(1,1);
goodspikes = [];

goodspikeind = (spikes(:,3) ~= -1) & (spikes(:,5) >= minabsvel);
%create a list of the non sharp-wave spikes [time traj linearloc] from
goodspikes = spikes(goodspikeind,:);
%make a cell array, where each cell contains data for one trajectory.
%inside each cell [binlocation occupancy spikecount firingrate]

trajdata = cell(1,trajnum);
goodlocationind = find((statevector ~= -1) & (absvel >= minabsvel));
goodlocations = [statematrix.time(goodlocationind) lindist(goodlocationind) statevector(goodlocationind)]; %the linear locations at all valid times


%     if isempty(goodspikes) %there were no useable spikes
%         goodspikes = [];
%         trajdata = [];
%     end

trajval = []; % trajdata has 4 possible trajectories for W-track
for i = 1:length(trajdata)  
    
    %get all the linear locations when the animal was on the ith
    %trajectory
    tmplinloc = goodlocations(find(goodlocations(:,3) == i),2); % lindist for current traj
    tmpspikes = goodspikes(find(goodspikes(:,3) == i),:); % Spikes for current traj
    if ~isempty(tmplinloc)
        findtrajnum = (i+rem(i,2))/2;
        minloc = 0;
        maxloc = trajlength(findtrajnum);
        
        %maxloc = max(tmplinloc);
        %minloc = min(tmplinloc);
        
        %numspatialbins = ceil((maxloc-minloc)/binsize));
        
        %binvector = zeros(numspatialbins,5);
        %locstep = (maxloc-minloc)/numspatialbins;
        %tmpvec = [minloc:locstep:maxloc]';
        tmpvec = [minloc:binsize:maxloc];
        binvector = zeros(length(tmpvec),5);
        %the 1st column of binvector is the location for that bin
        binvector(:,1) = tmpvec(1:end);
        %find which bins the linear locations and spikes fall into
        binind = lookup(tmplinloc,binvector(:,1));
        if ~isempty(tmpspikes)
            spikebinind = lookup(tmpspikes(:,4),binvector(:,1));
        else
            spikebinind = [];
        end
        %sum up the occupancy and spikes in each bin
        for j = 1:size(binvector,1)
            binvector(j,2) = sum(binind == j) * timestep;  %the second column is occupancy
            binvector(j,3) = sum(spikebinind == j); %3rd column is spike count
        end
        binvector(:,4) = binvector(:,3)./binvector(:,2); %4th column is occ-normalized firing rate
        nonfinite = find(~isfinite(binvector(:,4)));
        %binvector(nonfinite,4) = 0; %anything resulting from a divde-by-zero should be zero (there may be a better way to do this, because this causes the firing rates near zero-occupancy bins to be incorrectly lowered);
        %binvector(:,5) = gaussSmooth(binvector(:,4),2); %5th column is smoothed occ-normalized firing rate
        
        binvector(:,6) = gaussSmooth(binvector(:,2),2); %6th column is smoothed occupancy
        binvector(:,7) = gaussSmooth(binvector(:,3),2); %7th column is smoothed spike count
        binvector(:,5) = binvector(:,7)./binvector(:,6); %5th column is smoothed occ-normalized firing rate
        lowocc = find(binvector(:,6) < binsize*.1); % very low occupancy bins should be excluded
        
        
        %after the smoothing, turn the firing rate bins that had low occupancy to nan's.
        binvector(nonfinite,4) = nan;
        binvector(lowocc,5) = nan;
        
        trajdata{i} = binvector;
    end
    
end
warning('ON','MATLAB:divideByZero');


function  out = gaussSmooth(vector, binrange)

paddinglength = round(binrange*2.5);
padding = ones(paddinglength,1);

out = smoothvect([padding*vector(1) ; vector; padding*vector(end)],gaussian(binrange,binrange*5));
out = out(paddinglength+1:end-paddinglength);