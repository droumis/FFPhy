function data = createadaptdata(index, excludeperiods, linpos, spikes, timestep)
%function data = createadaptdata(index, excludeperiods, linpos, spikes, timestep)
%   returns the data structure for adaptive estimation of position given the
%   index element of the linearized position and spikes. The output times are
%   spaced at timestep. Invalid times given by excludeperiods are removed.
%
%   Note that this is genearlly run within the filtering framework.
%
%   data has three fields:
%	time 	a list of all timesteps
%	linpos	the animal's linearized distance from the start point for each
%		time
%	traj	the trajectory the animal is on at each time

time = linpos{index(1)}{index(2)}.statematrix.time;
data.time = time(1):timestep:time(end);

lp = linpos{index(1)}{index(2)}.statematrix.lindist;
data.linpos = interp1(time, lp, data.time, 'linear');

seginfo = linpos{index(1)}{index(2)}.segmentInfo;

traj = linpos{index(1)}{index(2)}.statematrix.traj;
data.traj = interp1(time, traj, data.time, 'nearest');

% note that the trajectory definitions from the statematrix use -1 for invalid
% elements and 1-n for the valid elements.  The adaptive filtering code expects
% the valid elements to start at 0, so we subtract 1 from all values > 0
tmp = find(data.traj > 0);
data.traj(tmp) = data.traj(tmp) - 1;
data.traj(find(isExcluded(data.time, excludeperiods))) = -1;
data.spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);

