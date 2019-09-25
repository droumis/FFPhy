function [h] = plotlinpos(linpos)
% function [h] = plotlinpos(linpos)
% Plots the linear distance from the start well for each of the trajectories,
% coloring each trajectory uniquely
%  
% linposstruct is an element of the linear position structure (e.g.
% linpos{2}{2})
%
state = linpos.statematrix;
traj = state.traj;
ntraj = length(unique(traj)) - 1;

refwell = state.referenceWell(1);

% plot all the invalid points in grey
invp = find(traj == 0);
if (isempty(invp))
    ntraj = ntraj+1;
else
    plot(state.time(invp), state.linearDistanceToWells(invp,  refwell), ...
	 '.', 'Color', [.5 .5 .5], 'MarkerSize', 6);
end

hold on
colors = jet(ntraj);

for i = 1:ntraj
    p = find(traj == i);
    plot(state.time(p), state.linearDistanceToWells(p,  refwell), '.', ...
	 'Color', colors(i,:), 'MarkerSize', 6);
end
h = gcf;
