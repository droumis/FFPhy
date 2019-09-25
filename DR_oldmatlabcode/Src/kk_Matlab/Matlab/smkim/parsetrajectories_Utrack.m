function trajectories = parsetrajectories_lineartrack(linearizedpos,min_trip_duration,min_excursion_distance)
% THIS FUNCTION CONTAINS HARD-CODED CONSTANTS!
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

% define functions for determining whether the linearized position is in 
% a food well zone
zones(1).descript = 'right food well';
% within 10cm of the end of the right arm
zones(1).test = @(linpos) (linpos(:,2)==1) & (linpos(:,3)>140);
zones(2).descript = 'left food well';
% within 10cm of the end of the left arm
zones(2).test = @(linpos) (linpos(:,2)==3) & (linpos(:,3)>140);

% integer values indicate whether rat is within either food well zone;
% zero values indicate that the rat is away from the food wells
currentstate = 1 * zones(1).test(linearizedpos.data(:,2:3)) + ...
    2 * zones(2).test(linearizedpos.data(:,2:3)) ];

% tmpidx is a temporary buffer of indices into the posdata
tmpidx = [];

% keep a count of the number of trajectories
currentcount = 0;


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
                tstart = ts2str(smoothedpos.data(tmpidx(1),1));
                tend = ts2str(smoothedpos.data(tmpidx(end),1));
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
        tstart = ts2str(smoothedpos.data(tmpidx(1),1));
        tend = ts2str(smoothedpos.data(tmpidx(end),1));
        traj.smoothedpos_snippets(currentcount).tstart = tstart;
        traj.smoothedpos_snippets(currentcount).tend = tend;
    else
        % otherwise, remove the final (false, abortive) element of startzone
        traj.startzone(currentcount) = [];
    end
elseif numel(traj.endzone) < (currentcount-1)
    error('bug in this program: endzone array can not be more than one element shorter than startzone array')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now score the trajectories as correct or error - this is really easy on the linear track!
% keep track of the last-visited side arm
for i = 1:numel(traj.startzone)
    if ~isnan(traj.startzone(i)) & ~isnan(traj.endzone(i))
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
end


