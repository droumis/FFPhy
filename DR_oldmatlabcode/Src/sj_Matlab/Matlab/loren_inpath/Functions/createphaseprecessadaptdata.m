function [data] = createphaseprecessadaptdata(sind, tind, excludeperiods, linpos, spikes, theta, timestep)
%function data = createphaseprecessadaptdata(sind, tind, excludeperiods, linpos, spikes,
%   returns the data structure for phase precession adaptive estimation given the
%   sind element of the linearized position and spikes. The output times are
%   spaced at timestep. Invalid times given by excludeperiods are removed.
%
%   Note that this is genearlly run within the filtering framework.
%
%   data has three fields:
%	time 	a list of all timesteps
%	linpos	the animal's linearized distance from the start point for each
%		time
%	traj	the trajectory the animal is on at each time
%	phase	the phase of the theta rhythm at each time


time = linpos{sind(1)}{sind(2)}.statematrix.time;
data.time = time(1):timestep:time(end);

lp = linpos{sind(1)}{sind(2)}.statematrix.lindist;
data.linpos = interp1(time, lp, data.time, 'linear');

traj = linpos{sind(1)}{sind(2)}.statematrix.traj;
data.traj = interp1(time, traj, data.time, 'nearest');

phase = geteegdata(theta{tind(1)}{tind(2)}{tind(3)}, 'phase');
eegtimes = geteegtimes(theta{tind(1)}{tind(2)}{tind(3)});

data.phase = interp1(eegtimes, phase, data.time, 'linear');


% note that the trajectory definitions from the statematrix use -1 for invalid
% elements and 1-n for the valid elements.  The adaptive filtering code expects
% the valid elements to start at 0, so we subtract 1 from all values > 0
tmp = find(data.traj > 0);
data.traj(tmp) = data.traj(tmp) - 1;
data.traj(find(isExcluded(data.time, excludeperiods))) = -1;
data.spiketimes = spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.data(:,1);


