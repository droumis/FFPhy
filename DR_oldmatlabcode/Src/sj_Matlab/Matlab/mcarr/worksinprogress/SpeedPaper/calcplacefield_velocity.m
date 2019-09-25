% Try pulling out each pass for each trajectory and assigning determining
% how fast it it.

%Load variables
load '/data13/mcarr/Eig/Eigspikes01.mat'
load '/data13/mcarr/Eig/Eigcellinfo.mat'
load '/data13/mcarr/Eig/Eiglinpos01.mat'
load '/data13/mcarr/Eig/Eigpos01.mat'
load '/data13/mcarr/Eig/Eigripples01.mat'

index = [1 2];
cells = evaluatefilter(cellinfo{index(1)}{index(2)},'isequal($area,''CA1'') & $peakrate>3');

%Define options
cmperbin = 5;
velocity_cutoff = 5;

linpos = linpos{index(1)}{index(2)}.statematrix;
p = linpos.lindist;
time = linpos.time;
timestep = median(diff(time));
spikes = spikes{index(1)}{index(2)};
vel = pos{index(1)}{index(2)}.data(lookup(time,pos{index(1)}{index(2)}.data(:,1)),8);

%Apply excludetimes
riptimes = getripples([index(1) index(2)], ripples, cellinfo, 'cellfilter', 'isequal($area,''CA1'')','minstd',3);
linpos.traj{6}(find(isExcluded(time, riptimes(:,[1 2])))) = -1;

traj = cell(4,1);
for t = 1:4
    %Define times for each traj
    traj{t}.times = linpos.traj{6}==t;

    %Exclude rippletimes
    
    % Define each trajectories distance
    traj{t}.bin = floor(min(p(traj{t}.times))):cmperbin:max(p(traj{t}.times));
    
    % Compute occupancy and speed
    traj{t}.occupancy = hist(p(traj{t}.times),traj{t}.bin).*timestep;
    traj{t}.speed = vel(traj{t}.times);
    
    % Determine fast and slow times and compute occupancy
    traj{t}.slow_occupancy = hist(p(traj{t}.times & vel<velocity_cutoff),traj{t}.bin).*timestep;
    traj{t}.fast_occupancy = hist(p(traj{t}.times & vel>velocity_cutoff),traj{t}.bin).*timestep;
    
    % Initialize each cell for each trajectory
    traj{t}.spikes = cell(size(cells,1),1);
    traj{t}.firing_slow = cell(size(cells,1),1);
    traj{t}.firing_fast = cell(size(cells,1),1);
    traj{t}.firing_total = cell(size(cells,1),1);
end

% Go through each cell for each pass and compute occupancy normalize firing
% rate
for c = 1:size(cells,1)
    for t = 1:4
        %Compute firing rate for all times
        validspikes = zeros(size(traj{t}.times));
        validspikes(spikes{cells(c,1)}{cells(c,2)}.data(:,7)) = 1;
        traj{t}.spikes{c} = hist(p(logical(validspikes)),traj{t}.bin);

        %Compute smoothed occupancy
        tmpocc = gaussSmooth_tmp(traj{t}.occupancy',2);
        %Compute smoothed spikecount
        tmpspike = gaussSmooth_tmp(traj{t}.spikes{c}',2);
        %Compute smoothed occupancy normalized spikecount
        traj{t}.firing_total{c} = tmpspike./tmpocc;
        traj{t}.firing_total{c}(traj{t}.occupancy < 0.1) = nan;
        
        %Compute firing rate for slow times
        slow_spikes = validspikes & vel < velocity_cutoff;
        tmpspike = hist(p(slow_spikes),traj{t}.bin);

        %Compute smoothed occupancy
        tmpocc = gaussSmooth_tmp(traj{t}.slow_occupancy',2);
        %Compute smoothed spikecount
        tmpspike = gaussSmooth_tmp(tmpspike',2);
        %Compute smoothed occupancy normalized spikecount
        traj{t}.firing_slow{c} = tmpspike./tmpocc;
        traj{t}.firing_slow{c}(traj{t}.slow_occupancy < 0.1) = nan;

        %Compute firing rate for fast times
        fast_spikes = validspikes & vel > velocity_cutoff;
        tmpspike = hist(p(fast_spikes),traj{t}.bin);

        %Compute smoothed occupancy
        tmpocc = gaussSmooth_tmp(traj{t}.fast_occupancy',2);
        %Compute smoothed spikecount
        tmpspike = gaussSmooth_tmp(tmpspike',2);
        %Compute smoothed occupancy normalized spikecount
        traj{t}.firing_fast{c} = tmpspike./tmpocc;
        traj{t}.firing_fast{c}(traj{t}.fast_occupancy < 0.1) = nan;
    end
end

% Go through and plot the firing rate map for the fast and slow times for
% each cell

figure
for c = 1:size(cells,1)
    for t = 1:4
        if max(traj{t}.firing_total{c}) > 3
            plot(traj{t}.bin, traj{t}.firing_total{c},'k',...
                traj{t}.bin, traj{t}.firing_slow{c}, 'r',...
                traj{t}.bin, traj{t}.firing_fast{c}, 'g')
            title(['Trajectory ', num2str(t),' Cell #', num2str([cells(c,1) cells(c,2)])])
            pause(5)
        end
    end
end
