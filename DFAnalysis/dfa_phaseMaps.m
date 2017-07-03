function [out] = dfa_phaseMaps(index, excludeperiods, spikes,linpos, pos, task, varargin)

%
% returns the theta phase angle pos, linpos, vel for each spike of each neuron
% organizes the output into the different trajectories.
%
%
%spikes - the 'spikes' cell array for the day you are analyzing
%linpos - the output of LINEARDAYPROCESS for the day you are analyzing.
%pos - the output of nspike_fixpos
%index - [day epoch tetrode cell]
%statevector - the outputs of GETBEHAVESTATE. This is a vector with the traj
%              number for each position time (1 based). -1 values signify
%              invalid times and are not used.


% parse the options
animals = [];
appendindex = 1;
binsize = 1; % cm squarec
stdocc = 3;
stdspike = 3;
threshocc = 0.002; % Threshold occupancy in seconds
if (~isempty(varargin))
    assign(varargin{:});
end

warning('OFF','MATLAB:divideByZero');
out.trajspikes = [];
out.allseqX = [];
out.allseqY = [];
out.allseqX1 = [];
out.allseqY1 = [];
out.wellCoords = [];
out.index = index;

try spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
catch
    %     output.spikes = [];
    %     out = output;
    return
end

spikesData = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
posdata = pos{index(1)}{index(2)}.data;
posfields = pos{index(1)}{index(2)}.fields;
taskEnv = task{index(1)}{index(2)}.environment;

timestring = 'time';
timecol = find(cell2mat(cellfun(@(x) strcmp(x,timestring), strsplit(posfields, ' '), 'UniformOutput', false)));
postime = posdata(:,timecol);

if strcmp(taskEnv, 'wtrack')
    statematrix = linpos{index(1)}{index(2)}.statematrix;
elseif strcmp(taskEnv, 'openfield')
    statematrix.time = postime;
    statematrix.traj = ones(length(postime),1);
end

statevector = statematrix.traj;

%DR added 4/22/14... front padding the posdata vect to match statevec
if length(statevector(:,1)) ~= length(posdata(:,1));
    posdata = [zeros((length(statevector(:,1))-length(posdata(:,1))), length(posdata(1,:))); posdata];
end

% Use Exclude Periods for TimeFilter version in addition to statevector=-1
statevector(find(isExcluded(posdata(:,1), excludeperiods))) = -1; % Based on exclude time,

posindexfield = knnsearch(postime, spikesData(:,1));

xstring = 'x-sm';
xcol = find(cell2mat(cellfun(@(x) strcmp(x,xstring), strsplit(posfields, ' '), 'UniformOutput', false)));
posxdata = posdata(posindexfield,xcol);
ystring = 'y-sm';
ycol = find(cell2mat(cellfun(@(x) strcmp(x,ystring), strsplit(posfields, ' '), 'UniformOutput', false)));
posydata = posdata(posindexfield,ycol);

if ~isempty(spikesData)
    spikesData = [spikesData(:,1) posxdata posydata posindexfield];
    spikesData(:,5) = statevector(posindexfield); %add the traj state for each spike
    %posindexfield = isdatafield(spikesfields,'posindex');
    
else
    %     spikesData = [0 0 -1];
    return
end

trajnum = max(statevector);
timestep = statematrix.time(2,1) - statematrix.time(1,1);
goodspikes = [];
goodspikeind = (spikesData(:,5) ~= -1);
%create a list of the non sharp-wave spikes [time traj linearloc]
goodspikes = spikesData(goodspikeind,:);
%make a cell array, where each cell contains data for one trajectory.
%inside each cell [binlocation occupancy spikecount firingrate]

trajdata = cell(1,trajnum);
goodlocationind = (find(statevector ~= -1 ));
goodlocations = [statematrix.time(goodlocationind) posdata(goodlocationind,[2 3]) statevector(goodlocationind)]; %CHANGED the positionindex at all valid times

out.allseqX(:,1) = 1:(((ceil(max(goodlocations(:,2)))-floor(min(goodlocations(:,2))))/binsize)+1);
out.allseqY(:,1) = 1:(((ceil(max(goodlocations(:,3)))-floor(min(goodlocations(:,3))))/binsize)+1);
out.allseqX(:,2) = floor(min(goodlocations(:,2))):binsize:ceil(max(goodlocations(:,2)));
out.allseqY(:,2) = floor(min(goodlocations(:,3))):binsize:ceil(max(goodlocations(:,3)));
out.allseqX1(:,2) = 1:binsize:ceil(max(goodlocations(:,2)));
out.allseqY1(:,2) = 1:binsize:ceil(max(goodlocations(:,3)));

animInfo = animaldef(lower(animals{an}));
% tmptheta = loadeegstruct(animInfo{2},animInfo{1},'theta',index(1),index(2),index(3));
tmptheta = loadeegstruct(animInfo{2},animInfo{1},'theta',index(1),index(2),11);
% thetaeeg = tmptheta{index(1)}{index(2)}{index(3)};
thetaeeg = tmptheta{index(1)}{index(2)}{11};
thetaeegtime = geteegtimes(thetaeeg)';
%     PhaseFreq=eegfilt(eeg_phs',samprate,Pf1,Pf2); % this is just filtering
% thetaphasetimeseries = angle(hilbert(double(thetaeeg.data(:,1)))); % this is getting the phase time series


% Resampled LFP oscillation phase
% It's very important to perform all computations in double precision, because
% MATLAB's mod function accumulates floating-point errors
% Map values to the range [-pi,+pi] and then unwrap
% call INTERP1 on the unwrapped angles and backtransform to
% recover values in [-pi,+pi]
thetaphasetimeseries = angle(hilbert(double(thetaeeg.data(:,1))));
% thetaphasetimeseries = thetaeeg.data(:,2);
spike2eegIndex = knnsearch(thetaeegtime, spikesData(:,1));
spikesData(:,6) = thetaphasetimeseries(spike2eegIndex);

% x = double(thetaeegtime);
% xi = double(spikesData(:,1));
% z = double(thetaphasetimeseries);
% a = mod(z+pi, 2*pi)-pi;
% b = unwrap(a,pi,1);
% c = interp1(x,b,xi);
% d = mod(c+pi, 2*pi)-pi;
% spikesData(:,6) = d;

% thetaphasetimeseries = thetaeeg.data(:,2);
% EEGindexfield = knnsearch(thetaeegtime, spikesData(:,1));
% spikethetaphase = thetaphasetimeseries(EEGindexfield);
% spikesData(:,6) = spikethetaphase/10000;
if strcmp(taskEnv, 'wtrack')
    lindistspike = statematrix.lindist(posindexfield);
    spikesData(:,7) = lindistspike;
end

for itraj = 1:length(trajdata)
    itrajspikesind = find(spikesData(:,5) == itraj);
    itrajspikes = spikesData(itrajspikesind,:);
    %add position and stuff to this traj output
    out.trajspikes{itraj} = itrajspikes;
    if strcmp(taskEnv, 'wtrack')
        wellCoords = linpos{index(1)}{index(2)}.wellSegmentInfo.wellCoord;
        traj2wells = [1 1 2; 2 2 1; 3 1 3; 4 3 1]; % traj enterwellID exitwellID
        out.wellCoords{itraj} = wellCoords(traj2wells(itraj,2:3),:); % to do [enterX enterY ; exitX exitY]
    else
        out.wellCoords{itraj} = [];
    end
end
disp(sprintf('day %d ep %d nt %d cell %d',index(1),index(2),index(3),index(4)))

if appendindex == 1
    out.index = index;
end

warning('ON','MATLAB:divideByZero');


