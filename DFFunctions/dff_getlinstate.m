function out = dff_getlinstate(animaldir,animalprefix,epochs, includeStates, varargin)

% Timefilter
% Shantanu - Equivalent to sj_getbehavestate in timefilter format. Returns
% additional parameters

% out = getalllinstate(animaldir,animalprefix,epochs, includeStates,options)
% Produces a cell structure with the fields:
% traj=state, lindist = lindist
% Also returns time, headdir, linearvelocity
% Head direction and linear velocity are linearized values.
% EPOCHS - N by 2 matrix, columns are [day epoch]
% INCLUDESTATES is a number between 1 and 6.  A higher 
%number tells this program to include progressively more ambiguous behavior
%times into statevector.  (1 is least ambiguous). Here are the definitions
%for each number:
%
%1:  Include all times when the animal is in transit between the two wells
%    of a defined trajectory (either in the foreward or backward direction),  
%    and velocity is higher than minvelocity. This is defined as the times 
%    when the animal is in transit between the two different endpoints and is
%    on a track segment that is part of the trajectory. Also, head
%    direction and motiondir have to be in the same direction.
%2:  Include times not occurring during any of the defined well-to-well trajectories, 
%    and occur when the animal is not leaving and entering the same well,
%    and is in a track segment that only belongs to one of the defined linear
%    trajectories. Also, head dir and motiondir
%    have to be in the same direction.
%3:  Fill all times that occur when the animal is on a segment 
%    that belongs to only one trajectory (even if it is leaving and entering the same well).
%    Also, head dir and motiondir have to be in the same direction.
%4:  Include times that occur when the animal was on an amiguous segment (like the
%    home arm). Also, head dir and motiondir have to be in the same direction.
%5:  Include any remaining points occuring in the beginning and end of the
%    session. Also, head dir and motiondir have to be in the same direction.
%6:  Include all points by choosing the closest defined trajectory
%
%
%Options:
%   'headir': Turn on or off the criterion that headdir and motion need to
%    be in the same direction. (default 1)
%   'minlinvelocity' - the minimum linear speed the animal must be
%           moving for the behavior to be counted as valid (default 6cm/sec)
%
%

loaddays = unique(epochs(:,1));
linpos = loaddatastruct(animaldir, animalprefix, 'linpos', loaddays); % All days in one linpos
for i = 1:size(epochs,1)
    lineindex = [1:length(linpos{epochs(i,1)}{epochs(i,2)}.statematrix.referenceWell)]';
    matrixIndex = sub2ind(size(linpos{epochs(i,1)}{epochs(i,2)}.statematrix.segmentHeadDirection),lineindex,linpos{epochs(i,1)}{epochs(i,2)}.statematrix.referenceWell);
    out{epochs(i,1)}{epochs(i,2)}.time = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.time;
    out{epochs(i,1)}{epochs(i,2)}.headdir = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.segmentHeadDirection(matrixIndex);
    out{epochs(i,1)}{epochs(i,2)}.linearvel = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.linearVelocity(matrixIndex);
    
    % Check if linpos has been updated to include traj=state and lindist
    if (isfield(linpos{epochs(i,1)}{epochs(i,2)}.statematrix,'traj')) && (isfield(linpos{epochs(i,1)}{epochs(i,2)}.statematrix,'lindist'))
        out{epochs(i,1)}{epochs(i,2)}.state = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.traj;
        out{epochs(i,1)}{epochs(i,2)}.lindist = linpos{epochs(i,1)}{epochs(i,2)}.statematrix.lindist;       
    else
  
        % Only recalculate state if state file does not exist
        % Existing state file is created by sj_getbehavestate
        statefile = sprintf('%s/%sbehavestate%02d.mat', animaldir, animalprefix, epochs(i,1));
        if ~exist(statefile,'file');
            [out{epochs(i,1)}{epochs(i,2)}.state, out{epochs(i,1)}{epochs(i,2)}.lindist] = getlinstate_internal(linpos, epochs(i,1), epochs(i,2), includeStates, varargin{:});
            %[out{epochs(i,1)}{epochs(i,2)}.state, out{epochs(i,1)}{epochs(i,2)}.lindist] = sj_getbehavestate(linpos{epochs(i,1)}, epochs(i,1), epochs(i,2), includeStates, varargin{:});
        else
            load(statefile);
            out{epochs(i,1)}{epochs(i,2)}.state = behavestate{epochs(i,1)}{epochs(i,2)}.state;
            out{epochs(i,1)}{epochs(i,2)}.lindist = behavestate{epochs(i,1)}{epochs(i,2)}.lindist;
        end        
    end
end
%--------------------------------------------------------------------------

function [state, lindist] = getlinstate_internal(linpos, day, epoch, includeStates, varargin)


timerange = [-inf inf];
minvelocity = 5;
headdir = 1;

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'headdir'
            headdir = varargin{option+1};
        case 'minlinvelocity'
            minvelocity = varargin{option+1};
         case 'timerange'
            timerange = varargin{option+1};
        otherwise
            error(['Option', varargin{option},'not defined']);
    end
end

statematrix = linpos{day}{epoch}.statematrix;
segmenttable = linpos{day}{epoch}.segmenttable;
trajwells = linpos{day}{epoch}.trajwells;

traject = trajwells;
includeStates = 1:includeStates;

%if no reference well is assigned, it is not a valid point
%also, only points inside the designated time range are valid
validPoints = ((statematrix.referenceWell > 0) & (statematrix.time >= timerange(1)) & (statematrix.time <= timerange(2))); 

validPointsIndex = find(validPoints);
trajvector = ones(length(validPointsIndex),1)*-1;
statevector = ones(length(validPointsIndex),1)*-1;
referencewell = statematrix.referenceWell(validPointsIndex);
welltraj = statematrix.wellExitEnter(validPointsIndex,:); %exit and enter wells
segmentindex = statematrix.segmentIndex(validPointsIndex);
samewell = (welltraj(:,1) == welltraj(:,2)); %Which times are during trajectories between the same wells 
diffwell = (welltraj(:,1) ~= welltraj(:,2)); %Which times are during trajectories between different wells
trajcount = rowcount(statematrix.segmentIndex(validPointsIndex),segmenttable(:,3)); %how many of the linear trajectories does each segment belong to

matrixIndex = sub2ind(size(statematrix.linearDistanceToWells),validPointsIndex,referencewell);
lineardistance = statematrix.linearDistanceToWells(matrixIndex);

if headdir == 1  %if head and motion direction should be in same direction
    forewarddir = ((statematrix.segmentHeadDirection(matrixIndex) >= 0) & (statematrix.linearVelocity(matrixIndex) >= minvelocity)); %facing positve dir and moving positive dir
    backwarddir = ((statematrix.segmentHeadDirection(matrixIndex) <= 0) & (statematrix.linearVelocity(matrixIndex) <= -minvelocity));%facing negative dir and moving negative dir

elseif headdir == 0 %if head and motion direction do not need to be in same direction
    forewarddir = ( (statematrix.linearVelocity(matrixIndex) >= minvelocity)); %facing positve dir and moving positive dir
    backwarddir = ((statematrix.linearVelocity(matrixIndex) <= -minvelocity));%facing negative dir and moving negative dir
end

%LEVEL 1   
for trajnum = 1:size(traject,1)
    currtraj = traject(trajnum,:);
    
    %which times occur when the animal is intransit between the two wells
    %for the current trajectory
    inDefinedTraj(:,trajnum) = (((welltraj(:,1)==currtraj(1))&(welltraj(:,2)==currtraj(2))) | ((welltraj(:,1)==currtraj(2))&(welltraj(:,2)==currtraj(1))));        
    %which segments are in the current trajectory
    segmentsInTraj = segmenttable(find(segmenttable(:,1) == trajnum),3);
    
    %find all times when the animal was on the current trajectory (either
    %in the foreward or backward direction).  This is defined as the times 
    %when the animal is in transit between the two correct endpoints and is
    %on a track segment that is part of the trajectory.
    
    %first find the foreward direction times (both head and motion
    %direction in the positive direction, and velocity greater than
    %minvelocity)       
    findindex = find( inDefinedTraj(:,trajnum)& ismember(segmentindex,segmentsInTraj) & forewarddir);   
    %assign the current trajectory to those indeces
    trajvector(findindex) = 2*(trajnum)-1;
    if ismember(1,includeStates)
        statevector(findindex) = trajvector(findindex);
    end
    %then find the negative direction times
    findindex = find( inDefinedTraj(:,trajnum) & ismember(segmentindex,segmentsInTraj) & backwarddir);        
    trajvector(findindex) = 2*(trajnum);
    if ismember(1,includeStates)
        statevector(findindex) = trajvector(findindex);
    end
end

undefinedindex = (trajvector == -1); %still undefined times
inAnyDefinedTraj = (sum(inDefinedTraj,2) > 0); %These times occur when the animal is in transit between two wells of a defined trajectory

%LEVEL 2
%for the times that are still undefined and not occurring during any of the above well-to-well trajectories, 
%and occur when the animal is not leaving and entering the same well,
%and is in a track segment that only belongs to one of the defined linear
%trajectories, assign the proper trajectory 
findindex = find((undefinedindex) & (diffwell) & (trajcount == 1) & (forewarddir) & (~inAnyDefinedTraj)); 
trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1))-1;
if ismember(2,includeStates)
    statevector(findindex) = trajvector(findindex);
end
findindex = find((undefinedindex) & (diffwell) & (trajcount == 1) & (backwarddir) & (~inAnyDefinedTraj)); 
trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1));
if ismember(2,includeStates)
    statevector(findindex) = trajvector(findindex);
end

if ((nargout == 1) & (max(includeStates) <= 2))
    return
end

%LEVEL 3
%fill all times when the animal is on a segment that belongs to only one
%trajectory (even if it is leaving and entering the same well)
undefinedindex = (trajvector == -1); %still undefined times
findindex = find((undefinedindex) & (trajcount == 1) & (forewarddir)); 
trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1))-1;
if ismember(3,includeStates)
    statevector(findindex) = trajvector(findindex);
end
findindex = find((undefinedindex) & (trajcount == 1) & (backwarddir)); 
trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1));
if ismember(3,includeStates)
    statevector(findindex) = trajvector(findindex);
end

%LEVEL 4
%Fill in the undefined times when the animal was on an amiguous segment (like the
%home arm in a w-maze) 
undefinedindex = (trajvector == -1); %the undefined times 
%if facing positive direction, then assign the traj of the next traj in the
%future
findindex = find((undefinedindex) & (trajcount > 1) & (forewarddir)); 
trajvector = vectorfill(trajvector, -1, 1, findindex);
findex2 = findindex(find(~mod(trajvector(findindex),2)));
trajvector(findex2) = trajvector(findex2)-1;
if ismember(4,includeStates)
    statevector(findindex) = trajvector(findindex);
end
%otherwise assign the closest traj that happened in the past 
findindex = find((undefinedindex) & (trajcount > 1) & (backwarddir)); 
trajvector = vectorfill(trajvector, -1, -1, findindex);
findex2 = findindex(find(mod(trajvector(findindex),2)));
trajvector(findex2) = trajvector(findex2)+1;
if ismember(4,includeStates)
    statevector(findindex) = trajvector(findindex);
end

%LEVEL 5
% now we may have some undefined points in the beginning and end of the
% session- we assign these to the first trajectory
undefinedindex = (trajvector < 0); %still undefined times
trajvector(find(undefinedindex)) = -1;
statevector(find(undefinedindex)) = -1;
findindex = find((undefinedindex) & (forewarddir)); 
trajvector(findindex) = 1;
if ismember(5,includeStates)
    statevector(findindex) = trajvector(findindex);
end
findindex = find((undefinedindex) & (backwarddir)); 
trajvector(findindex) = 2;
if ismember(5,includeStates)
    statevector(findindex) = trajvector(findindex);
end

%LEVEL 6
% Fill in all remaining points, regardless of running velocity
undefinedindex = find(trajvector == -1); %still undefined times
if length(find(undefinedindex)) ~= length(trajvector)
    %trajvector = vectorfill(trajvector, -1, 0, undefinedindex);
    trajvector = vectorfill(trajvector, -1);
end
if ismember(6,includeStates)
    statevector = trajvector;
end


state = ones(length(statematrix.time),1)*-1;
lindist = ones(length(statematrix.time),1)*-1;
state(validPointsIndex) = statevector;
lindist(validPointsIndex) = lineardistance;
