% To plot linear position vs. occupancy normalized firing rates. First, use
% getbehavestate to find the linear trajectory for each time point (and
% determine invalid behavioral times). Then use calclinfields to calculate
% the occupancy normalized firing rates. Output is graph of linear position
% vs. occupancy normalized firing rate for outbound trajectory in red and
% return trajectory in blue.

% Define variables: day, epoch, tetrode, cell
day = 2
epoch = 6
tetrode = 2
cell = 1

% load variables
load (['m2linpos0', num2str(day)]);
load (['m2spikes0', num2str(day)]);

% Calculate occupancy normalized firing rates
[state, lindist] = getbehavestate (linpos, day, epoch, 2);
cell = calclinfields (spikes, state, lindist, linpos, [day epoch tetrode cell]);

% Plot figure
figure
plot (cell{1}(:,1) , cell{1}(:,5))
hold on
plot (cell{2}(:,1) , cell{2}(:,5), 'r')
xlabel('linear position (cm)')
ylabel('occupancy normalized firing rate')
title(['Day ', int2str(day), ', Epoch ', int2str(epoch), ', Tetrode ', int2str(tetrode)])