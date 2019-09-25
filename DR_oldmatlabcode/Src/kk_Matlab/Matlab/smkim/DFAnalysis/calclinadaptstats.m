function out = calclinadaptstats(index, excludeperiods, spikes, linpos, fout)

% out = calclinadaptstats(index, excludeperiods, spikes, linpos, fout)
%
% Calculates the statistics on the adaptive estimator.  The statistics and the
% locations where they should be computed are unique to each track, so you will
% need to edit this function for each different environment or analysis.
% 
% fout is the output of the filter from a previous run of calclinadapt
%

timestep = 0.002; % 2ms timestep
%timestep = fout.timestep;
data = createadaptdata(index, excludeperiods, linpos, spikes, timestep);

variables{1}.name = 'Moments'; % this produces four measures:  Area, mean, standard deviation and skewness
variables{1}.sizex = 1; %1 cm bins

selection.a{1}.traj= 0;           % trajectory to analyze (starting from 0)
selection.a{1}.linpos= [fout.x.cpx{1}(1) fout.x.cpx{1}(end)];    % entire trajectory
times = 1:fout.x.outputInterval:length(data.time);
selection.x{1}.time=  data.time(times); % times for analysis (sec)

out = calcStatistics(data, fout, variables, selection);
out.times = selection.x{1}.time;
out.index = index;



