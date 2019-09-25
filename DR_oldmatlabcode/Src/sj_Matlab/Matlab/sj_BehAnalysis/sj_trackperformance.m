function [trajbound, rewarded, time] = sj_trackperformance(index, excludetimes, linpos, correctorder, varargin)
% function [trajbound, rewarded, time] = trackperformance(index, excludetimes, linpos, correctorder, mindist)
%%assumes animal is performing an alternation task, with a middle well and
%2, alternating outside wells, though it does not have to be a W-shaped
%track.  does not track trajectories, as animal may wander between time
%leaves one well and enters another. 
%
% LINPOS is the output of linearizeposition for one day
% INDEX [day epoch]
% EXCLUDETIMES times to exclude from analysis.  This excludes trials whose
% start time falls in an exclude period
% CORRECTORDER is a vector of 3 well numbers, the middle number corresponds
% to the middle well, the first and third ones correspond to outside wells,
% for instance [2 1 3] would mean the middle well is well 1, the outside
% wells are 2 and 3.  these well numbers are determined in createtaskstruct
% varargin: MINDIST is the minimum length of an excursion from a well for the excursion
% to be counted as a trajectory 
%
% TRAJBOUND defined by where animal heading, 
%   -1 if animal heading to a well not included in correct order
%   1 if headed to middle well
%   0 if headed to outside well
% REWARDED, 0 or 1, 1 if rewarded at well heading to, -1 if not
% TIME time of traj start


mindist=[];
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'mindist'
            mindist = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

day = index(1);
epoch = index(2);
          
midwell=correctorder(2); %middle well
owell1 =correctorder(1); %1st outside well
owell2 =correctorder(3); %2nd outside well

wellExitEnter = linpos{day}{epoch}.statematrix.wellExitEnter;
startwell = wellExitEnter(:,1); %well leaving from
endwell = wellExitEnter(:,2); %well heading to
wellchange = [1; (find((diff(startwell) | diff(endwell))) + 1)];

% Shantanu - I already have turn-around errors by distance implements in linpos. 
% Can put in here again, but use proportion of length of track

% we need to get rid of cases where the trajectory indicates that the animal
% has started and ended at the same well but where the excursion distance is
% less than mindist

if ~isempty(mindist)
    samewell = find(wellExitEnter(wellchange,1) == wellExitEnter(wellchange,2));
    invalid = [];
    for s = 1:length(samewell)
        % check to find the maximum excursion distance for this set of points
        if (samewell(s) < length(wellchange))
            tind = wellchange(samewell(s)):wellchange(samewell(s)+1);
        else
            % this is the last trajectory, so we go to the end of the data
            tind = wellchange(samewell(s)):size(wellExitEnter,1);
        end
        % get the distances to the start well
        d = linpos{day}{epoch}.statematrix.linearDistanceToWells(tind, ...
            wellExitEnter(wellchange((samewell(s))),1));
        if (max(d) < mindist)
            invalid(s) = 1;
        else
            invalid(s) = 0;
        end
    end
    % get rid of the invalid trajectories
    wellchange(samewell(find(invalid))) = -1;
    wellchange = wellchange(find(wellchange ~= -1));
end

wellSequence = endwell(wellchange); %seq of end wells visited
correct = 0; %tracks number correct traj
rewarded = -5*ones(length(wellSequence),1); %tracks if each traj correct (1) or incorrect(0)
trajbound = rewarded;
total = 0;%tracks total number traj counted
lastoutwell = []; %tracks last outward bound well
time = linpos{day}{epoch}.statematrix.time(wellchange);
% we also need a list of end times.  This needs to be checked to make sure it
% is correct for the case of an incomplete last trajectory
endtime = [time(2:end) ; linpos{day}{epoch}.statematrix.time(end)];

%makes a vector: correctSeq of same length as wellSeq that is 0 if traj is
%incorrect, 1 if traj is incorrect
i = 1; %first traj
if ( ((wellSequence(i) == owell1) || (wellSequence(i) == owell2)) || (wellSequence(i) == midwell) ) 
    correct = correct+1;
    tmpR = 1;
    total = total+1;
    tmpT = 1;
    if ( ((wellSequence(i) == owell1) || (wellSequence(i) == owell2)) )
        tmpT= 0;
    end
else
    total = total+1;
    tmpR = 0;
    tmpT = -1;
end
rewarded(1) = tmpR;
trajbound(1) = tmpT;

for i = 2:length(wellSequence)
    tmpR = -5; tmpT = -5;
    if ( ((wellSequence(i) == owell1) || (wellSequence(i) == owell2)) & (wellSequence(i-1) == midwell) ) %is an outbound traj
        tmpT= 0;
        if i==2
            correct = correct+1;
            tmpR = 1;
            total = total+1;
        elseif ( (wellSequence(i-2) ~= wellSequence(i)) ) %if last outward bound well is different from current outward bound well    
            correct = correct+1;
            tmpR = 1;
            total = total+1;
        else
            total = total+1;
            tmpR = 0;
        end
    elseif ( ((wellSequence(i) == owell1) || (wellSequence(i) == owell2)) & (wellSequence(i-1) ~= midwell) )
        tmpR = 0;
        tmpT = 1;
    elseif ((wellSequence(i) == midwell) & (wellSequence(i-1) == midwell))
	%error
        tmpT = 0;
        tmpR = 0;
    elseif ( (wellSequence(i) == midwell) ) %inbound traj
        tmpT = 1;
        tmpR = 1;
    else 
        total = total+1;
        tmpR = 0;
        tmpT = -1;
    end

    if ( ((wellSequence(i) ==owell1) || (wellSequence(i) == owell2)) )
        lastoutwell = wellSequence(i); %set last outward well to current outwell for next trial
    end
    rewarded(i) = tmpR;
    trajbound(i) = tmpT;
end

%apply excludetimes
include = find(~isExcluded(time, excludetimes));
time = [time(include) endtime(include)];
rewarded = rewarded(include);
trajbound = trajbound(include);

