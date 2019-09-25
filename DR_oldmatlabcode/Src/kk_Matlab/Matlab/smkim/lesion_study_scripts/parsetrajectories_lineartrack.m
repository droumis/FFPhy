function traj = parsetrajectories_lineartrack(smoothedpos,min_trip_duration,min_excursion_distance)
% 
%   traj = parsetrajectories_lineartrack(smoothedpos,min_trip_duration,min_excursion_distance)
%
%       .descript
%       .subject
%       .day
%       .session
%       .environment
%       .tstart
%       .tend
%       .registeredframe
%       .template
%       .params
%       .min_trip_duration
%       .min_excursion_distance
%       .startzone(i)
%       .endzone(i)
%       .correct(i)
%       .smoothedpos_snippets(i)
%           .fields
%           .data
%           .tstart
%           .tend
%

% inherit metdata from smoothedpos
traj = smoothedpos;
traj.descript = 'trajectories on linear track';
traj.min_trip_duration = min_trip_duration;
traj.min_excursion_distance = min_excursion_distance;

% now we determine, at every position sample, whether the animal was inside 
% each food-well zone defined in the track template.
% each column of zone_occupancies is a list of flags; zone_occupancies(j,i) 
% indicates whether the animal was within zones(i) at pos sample (j)
zones = traj.template.zones;
zone_occupancies = zeros(size(smoothedpos.data,1),length(zones));
for i = 1:length(zones)
    zone_occupancies(:,i) = inpolygon(smoothedpos.data(:,2), ...
        smoothedpos.data(:,3), ...
        zones(i).vertices(:,1),zones(i).vertices(:,2));
end
% collapse these into a single vector listing the curent zone
currentzone = nan(size(zone_occupancies,1),1);
for j = 1:size(zone_occupancies,1)
    if nnz(zone_occupancies(j,:)) > 1
        error('you defined the zones so that they overlap!')
    end
    if any(zone_occupancies(j,:));
        currentzone(j) = find(zone_occupancies(j,:));
    else
        % implied NaN from initialization of currentzone vector
        % currentzone(j) = NaN;
    end
end
% tmpidx is a temporary buffer of indices into the posdata
tmpidx = [];

% initialize traj depending on where the rat starts
currentcount = 1;
lastzone = currentzone(1);
if ~isnan(lastzone) % if rat begins in a defined zone, continue look ahead 
    % until rat leaves this zone
    %pass
else % if rat begins not within a zone, pick up pos samples immediately
    traj.startzone(currentcount) = NaN;
    tmpidx = 1; % grab the first position sample
end

for j = 2:numel(currentzone)
    if ~isnan(lastzone) % if the rat was within a real zone in last pos sample
        if isnan(currentzone(j)) % if the rat leaves the zone
            traj.startzone(currentcount) = lastzone; % remember the last zone
            tmpidx = []; tmpidx = j; % clear tmpidx and grab a new start index
        elseif currentzone(j) == lastzone % if the rat stays in the zone
            % pass
        else
            error('rat jumped directly from one zone to another!');
        end
    else % if the rat was not within a defined zone at the last position sample
        if ~isnan(currentzone(j)) % if rat enters a new zone (i.e. approaches a food well)
            % check that the trajectory, starting from the last food well going to
            % this new food well, satisfies the min_trip_duration and 
            % min_excursion_distance criteria
            if ~isempty(tmpidx) % if we have posindices in our buffer
                % look up start and end times corresponding to tmpidx
                visit_duration = smoothedpos.data(tmpidx(end),1) - smoothedpos.data(tmpidx(1),1);
                excursion_dist = max(hypot( ...
                    smoothedpos.data(tmpidx(1),2) - smoothedpos.data(tmpidx(:),2), ...
                    smoothedpos.data(tmpidx(1),3) - smoothedpos.data(tmpidx(:),3) ));
            else
                visit_duration = 0;
                excursion_dist = 0;
            end
            % if the current excursion fails the test, overwrite the (currentcount) entry
            if (visit_duration < min_trip_duration) | ...
               (excursion_dist < min_excursion_distance)
                % reset the startzone
                traj.startzone(currentcount) = currentzone(j);
                tmpidx = []; % clear posindices buffer
            else % otherwise, cap the end of this trajectory
                traj.endzone(currentcount) = currentzone(j);
                traj.smoothedpos_snippets(currentcount).data = smoothedpos.data(tmpidx,:);
                traj.smoothedpos_snippets(currentcount).fields = smoothedpos.fields;
                tstart = timetrans(smoothedpos.data(tmpidx(1),1),1,1); tstart = tstart{1};
                tend = timetrans(smoothedpos.data(tmpidx(end),1),1,1); tend = tend{1};
                traj.smoothedpos_snippets(currentcount).tstart = tstart;
                traj.smoothedpos_snippets(currentcount).tend = tend;
                % advance to the next element of the currentcount struct array
                currentcount = currentcount + 1;
            end
        else % if the rat continues without a defined zone
            tmpidx(end+1) = j; % extend the posindices buffer
        end
    end 
    % remember the currentzone for the next iteration
    lastzone = currentzone(j);
end
currentcount = numel(traj.startzone);
% what do we do if final endzone is left undefined?
if numel(traj.endzone) == (currentcount-1)
    if ~isempty(tmpidx)
        % if we we have a non-empty tmpidx buffer, cap the final trajectory
        traj.endzone(currentcount) = NaN;
        traj.smoothedpos_snippets(currentcount).data = smoothedpos.data(tmpidx,:);
        traj.smoothedpos_snippets(currentcount).fields = smoothedpos.fields;
        tstart = timetrans(smoothedpos.data(tmpidx(1),1),1,1); tstart = tstart{1};
        tend = timetrans(smoothedpos.data(tmpidx(end),1),1,1); tend = tend{1};
        traj.smoothedpos_snippets(currentcount).tstart = tstart;
        traj.smoothedpos_snippets(currentcount).tend = tend;
    else
        % otherwise, remove the final (false, abortive) element of startzone
        traj.startzone(currentcount) = [];
    end
elseif numel(traj.endzone) < (currentcount-1)
    error('bug in this program: endzone array can not be more than one element shorter than startzone array')
end

%compute running speed over each journey
for i = 1:numel(traj.endzone)
    traj.speed(i) = mean(hypot(traj.smoothedpos_snippets(i).data(:,4), ...
        traj.smoothedpos_snippets(i).data(:,5)));
end

% compute duration of food-well visit at the end of every journey
for i = 1:numel(traj.endzone)
    if ~isnan(traj.startzone(i)) & ~isnan(traj.endzone(i)) & (i < numel(traj.endzone)) & ...
        (traj.smoothedpos_snippets(i).data(end,1) < traj.data(end,1))
        traj.visitduration(i) = ...
            traj.smoothedpos_snippets(i+1).data(1,1) - traj.smoothedpos_snippets(i).data(end,1);
    else
        traj.visitduration(i) = NaN;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now score the trajectories as correct or error - this is really easy on the linear track!
last_visited = traj.startzone(1);
for i = 1:numel(traj.startzone)
    if ~isnan(traj.startzone(i)) & ~isnan(traj.endzone(i)) & ~isnan(last_visited)
        % make sure that both startzone and endzone are defined
        if traj.endzone(i) == traj.startzone(i)
            % perseverative re-visit
            traj.correct(i) = 0;
        else
            % end-to-end running is correct
            traj.correct(i) = 1;
        end
    else
        traj.correct(i) = NaN;
    end
    last_visited = traj.endzone(i);
end


