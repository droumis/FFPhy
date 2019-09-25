function out = DFTFsj_getvelpos(animaldir,animalprefix,epochs, varargin)

% Timefilter
% Shantanu - Get Speed and position from animal from pos file

% out = DFsj_getvelpos(animaldir,animalprefix,epochs, options)

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
pos = loaddatastruct(animaldir, animalprefix, 'pos', loaddays); % All days in one linpos

for i = 1:size(epochs,1)
    
    % Position Data with Fields 'time x y dir vel x-sm y-sm dir-sm vel-sm'
    posdata = pos{epochs(i,1)}{epochs(i,2)}.data;
    if size(posdata,2)>5 % already smoothed position and filtered velocity
        abspos = posdata(:,6:7);
        absvel = abs(posdata(:,9)); % Should already be absolute value
    end
    % Using original position and velocity
    abspos = posdata(:,2:3);
    absvel = abs(posdata(:,5));
    
    out{epochs(i,1)}{epochs(i,2)}.absvel = absvel;
    out{epochs(i,1)}{epochs(i,2)}.abspos = abspos;
    out{epochs(i,1)}{epochs(i,2)}.time = posdata(:,1);
end

%--------------------------------------------------------------------------
