function linpos = linearizepos_Wtrack_old(smoothedpos,varargin)
%
%   linpos = linearize_Wtrack_old(smoothedpos,Wtrackdefinition)
%
% This function attempts an automated linearization over as much 
% of the data as it can manage. It then provides a gui for the 
% user to pick from among candidate linearizations at position 
% samples where the algorithm can not determine an unambiguous 
% linearization. 
%
% Wtrackdefinition is a standard Wtrack struct, with fields 
% .vertices and .posrange. If this argument is not included in 
% the function call, then the function defineWtrack will be 
% called so that the user can provide one.
%
% The output is a linpos.
%

% Definition of function constants
% Nominal frame rate in fps
FRAMERATE = 30;
% compare empirically observed timestamps diffs to the expected
if any(abs(FRAMERATE*diff(smoothedpos.data(:,1)) - 1) > 0.1)
    error('position samples deviate significantly from expected framerate')
end

% what is the maximum possible jump in linpos over a single frame? 
MAXJUMP = 200/FRAMERATE; % cm/frame

% proximity parameter for deciding how close is "close enough"
% this should be a small value (smaller will be more
% conservative, less error-prone)
PROXIMITY = 2; % 2 cm

% Error checking on the posparams input
if any(isnan(smoothedpos.data(:,4:5)))
    error('posparams data contains NaN values; you need to first interpolate over these');
elseif any(diff(smoothedpos.data(:,1)) =< 0)
    error('posparams data contains out-of-order or duplicate timestamps');
elseif any(hypot( ...
       diff(smoothedpos.data(:,4)),diff(smoothedpos.data(:,5))) > MAXJUMP)
    error('posparams data must not contain any large jumps in position');
end

% -----------------------------------------------------------------
% Construct a large collection of vectors for partitioning the trajectory
% into zones with different linearization schemes
v_a = Wtrack.vertices(1,:);
v_b = Wtrack.vertices(2,:);
v_c = Wtrack.vertices(3,:);
v_d = Wtrack.vertices(4,:);
v_e = Wtrack.vertices(5,:);
v_f = Wtrack.vertices(6,:);

% figure out the orientation of the W track with respect to the coordinate frame
% and define the bounding corners accordingly
posrange = [ ...
    min(smoothedpos.data(:,4))-20, ...
    max(smoothedpos.data(:,4))+20, ...
    min(smoothedpos.data(:,5))-20, ...
    max(smoothedpos.data(:,5))+20 ];
[distance, nearest_corner] = min( hypot( ...
    v_a(1) - [posrange(1); posrange(2); ...
    posrange(2); posrange(1)], ...
    v_a(2) - [posrange(3); posrange(3); ...
    posrange(4); posrange(4)]) );
switch nearest_corner
    case 1
        v_RightBottom = [ posrange(1) posrange(3) ];
        v_RightTop = [ posrange(2) posrange(3) ];
        v_LeftTop = [ posrange(2) posrange(4) ];
        v_LeftBottom = [ posrange(1) posrange(4) ];
    case 2
        v_RightBottom = [ posrange(2) posrange(3) ];
        v_RightTop = [ posrange(2) posrange(4) ];
        v_LeftTop = [ posrange(1) posrange(4) ];
        v_LeftBottom = [ posrange(1) posrange(3) ];
    case 3
        v_RightBottom = [ posrange(2) posrange(4) ];
        v_RightTop = [ posrange(1) posrange(4) ];
        v_LeftTop = [ posrange(1) posrange(3) ];
        v_LeftBottom = [ posrange(2) posrange(3) ];
    case 4
        v_RightBottom = [ posrange(1) posrange(4) ];
        v_RightTop = [ posrange(1) posrange(3) ];
        v_LeftTop = [ posrange(2) posrange(3) ];
        v_LeftBottom = [ posrange(2) posrange(4) ];
    otherwise
        error('There is something wrong with the Wtrack definition');
end

% define vectors in the directions along the track segments; note that
% this also provides information on the lengths of the track segments
v_1 = v_a - v_d;
v_2 = v_d - v_e;
v_3 = v_b - v_e;
v_4 = v_f - v_e;
v_5 = v_c - v_f;
% if segments 2 and 4 differ in length by more than 5cm, raise an error
if norm(v_2) - norm(v_4) > 5
    error('your track definition is such that track segment lengths are grossly unequal!');
end

% Linpos coordinates convention:
%
% segment 1: lindist is measured from vertex "d"
% segment 2: lindist is measured from vertex "e"
% segment 3: lindist is measured from vertex "e"
% segment 4: lindist is measured from vertex "e"
% segment 5: lindist is measured from vertex "f"
%
%          f--4--e--2--d     
%          |     |     |     
%          |     |     |     
%          5     3     1     
%          |     |     |     
%          |     |     |     
%          c     b     a     
%

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!! VERY IMPORTANT THAT segments BE DECLARED GLOBAL !!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
global segments;
segments = Wtrack.segments;

% Locate a point along segment 3 at which to draw a perpendicular
% for partitioning the space
v_OriginCenter = segments(3).origin + ... 
    1/4*(segments(2).length + segments(4).length) * ... 
    segments(3).direction;

% To find a vector which bisects the angle between two vectors: 
% (1) normalize each argument vector; (2) add the normalized
% unit vectors, paying attention to sign (direction)
v_SplitRight = v_1/norm(v_1) + v_3/norm(v_3);
v_SplitLeft = v_5/norm(v_5) + v_3/norm(v_3);

if norm( v_2/norm(v_2) + v_4/norm(v_4) ) == 0
    v_SplitCenter = perp(v_2);
else
    v_SplitCenter = v_2/norm(v_2) + v_4/norm(v_4);
end
if abs(dot(v_3,v_1)) == norm(v_3)*norm(v_1)
    v_OriginRight = ... % in case the arms are parallel
    0.5*v_d + 0.5*xlines(v_d,v_e,perp(v_1),v_3);
else
    v_OriginRight = xlines(v_d,v_e,v_1,v_3);
end
if abs(dot(v_3,v_5)) == norm(v_3)*norm(v_5)
    v_OriginLeft = ... % in case the arms are parallel
    0.5*v_f + 0.5*xlines(v_f,v_e,perp(v_5),v_3);
else
    v_OriginLeft = xlines(v_f,v_e,v_5,v_3);
end

% Define key "pivot" points for linearizing around the turns
v_PivotRight = xlines( ...
    v_OriginCenter,v_OriginRight, ...
    perp(v_3),v_SplitRight);

v_PivotLeft = xlines( ...
    v_OriginCenter,v_OriginLeft, ...
    perp(v_3),v_SplitLeft);

% Now that all of the key points and direction vectors are defined, 
% define vertices for the zones of interest which tile the space 
% (note that some zones overlap - this is by design). Depending on
% which zone(s) the trajectory enters, different linearization schemes
% will be applied. We also define for each zone the corresponding linear 
% segment of the track onto which the position is to be projected. Lastly,
% we define for each zone an anonymous function which computes the
% linearization from xypos ([x y]) to linpos ([segment distance])
zones(1).vertices = [ ...
    xlines(v_PivotRight,v_RightBottom,v_SplitRight,v_LeftBottom-v_RightBottom); ...
    v_RightBottom; ...
    xlines(v_PivotRight,v_RightBottom,perp(v_1),v_RightTop-v_RightBottom); ...
    v_PivotRight ];
zones(1).linearization = @(xypos) [1 dot(xypos-segments(1).origin,segments(1).direction)];

zones(2).vertices = [ ...
    v_PivotRight; ...
    xlines(v_PivotRight,v_RightBottom,perp(v_1),v_RightTop-v_RightBottom); ...
    xlines(v_PivotRight,v_RightTop,v_d-v_PivotRight,v_RightBottom-v_RightTop); ...
    v_RightTop; ...
    xlines(v_PivotRight,v_RightTop,v_d-v_PivotRight,v_LeftTop-v_RightTop) ];
zones(2).linearization = @(xypos) [1 norm( ...
    xlines(xypos,segments(1).origin,xypos-v_PivotRight,segments(1).direction) - ...
    segments(1).origin )];
     
zones(3).vertices = [ ...
    v_PivotRight; ...
    xlines(v_PivotRight,v_RightTop,v_d-v_PivotRight,v_RightBottom-v_RightTop); ...
    v_RightTop; ...
    xlines(v_PivotRight,v_RightTop,v_SplitRight,v_LeftTop-v_RightTop); ...
    xlines(v_PivotRight,v_RightTop,v_SplitRight,v_LeftTop-v_RightTop) ];
zones(3).linearization = @(xypos) [2 norm( ...
    xlines(xypos,segments(2).origin,xypos-v_PivotRight,segments(2).direction) - ...
    segments(2).origin )];

zones(4).vertices = [ ...
    v_PivotRight; ...
    xlines(v_PivotRight,v_RightTop,v_SplitRight,v_LeftTop-v_RightTop); ...
    xlines(v_PivotLeft,v_RightTop,v_SplitLeft,v_LeftTop-v_RightTop); ...
    xlines(v_PivotLeft,v_e,v_SplitLeft,v_e-v_PivotRight) ];
zones(4).linearization = @(xypos) [2 norm( ...
    xlines(xypos,segments(2).origin,xypos-v_PivotRight,segments(2).direction) - ...
    segments(2).origin )];
    
zones(5).vertices = [ ...
    v_PivotRight; ...
    xlines(v_PivotLeft,v_e,v_SplitLeft,v_e-v_PivotRight); ...
    v_PivotLeft ];
zones(5).linearization = @(xypos) [3 norm( ...
    xlines(xypos,segments(3).origin,xypos-v_PivotRight,segments(3).direction) - ...
    segments(3).origin )];

zones(6).vertices = [ ...
    xlines(v_PivotLeft,v_LeftBottom,v_SplitLeft,v_RightBottom-v_LeftBottom); ...
    xlines(v_PivotRight,v_RightBottom,v_SplitRight,v_RightBottom-v_LeftBottom); ...
    v_PivotRight; ...
    v_PivotLeft ];
zones(6).linearization = @(xypos) [3 dot(xypos-segments(3).origin,segments(3).direction)];

zones(7).vertices = [ ...
    v_PivotLeft; ...
    v_PivotRight; ...
    xlines(v_PivotRight,v_e,v_SplitRight,v_e-v_PivotLeft) ];
zones(7).linearization = @(xypos) [3 norm( ...
    xlines(xypos,segments(3).origin,xypos-v_PivotLeft,segments(3).direction) - ...
    segments(3).origin )];

zones(8).vertices = [ ...
    v_PivotLeft; ...
    xlines(v_PivotRight,v_e,v_SplitRight,v_e-v_PivotLeft); ...
    xlines(v_PivotRight,v_RightTop,v_SplitRight,v_LeftTop-v_RightTop); ...
    xlines(v_PivotLeft,v_RightTop,v_SplitLeft,v_LeftTop-v_RightTop) ];
zones(8).linearization = @(xypos) [4 norm( ...
    xlines(xypos,segments(4).origin,xypos-v_PivotLeft,segments(4).direction) - ...
    segments(4).origin )];

zones(9).vertices = [ ...
    v_PivotLeft; ...
    xlines(v_PivotLeft,v_LeftTop,v_SplitLeft,v_LeftTop-v_RightTop); ...
    xlines(v_PivotLeft,v_LeftTop,v_f-v_PivotLeft,v_LeftTop-v_RightTop); ...
    v_LeftTop; ...
    xlines(v_PivotLeft,v_LeftTop,v_f-v_PivotLeft,v_LeftBottom-v_LeftTop) ];
zones(9).linearization = @(xypos) [4 norm( ...
    xlines(xypos,segments(4).origin,xypos-v_PivotLeft,segments(4).direction) - ...
    segments(4).origin )];

zones(10).vertices = [ ...
    v_PivotLeft; ...
    xlines(v_PivotLeft,v_LeftTop,v_f-v_PivotLeft,v_LeftTop-v_LeftBottom); ...
    xlines(v_PivotLeft,v_LeftTop,v_f-v_PivotLeft,v_LeftBottom-v_RightTop); ...
    v_LeftTop; ...
    xlines(v_PivotLeft,v_LeftBottom,perp(v_5),v_LeftTop-v_LeftBottom) ];
zones(10).linearization = @(xypos) [5 norm( ...
    xlines(xypos,segments(5).origin,xypos-v_PivotLeft,segments(5).direction) - ...
    segments(5).origin )];

zones(11).vertices = [ ...
    v_LeftBottom; ...
    xlines(v_PivotLeft,v_LeftBottom,v_SplitLeft,v_RightBottom-v_LeftBottom); ...
    v_PivotLeft; ...
    xlines(v_PivotLeft,v_LeftBottom,perp(v_5),v_LeftTop-v_LeftBottom) ];
zones(11).linearization = @(xypos) [5 dot(xypos-segments(5).origin,segments(5).direction)];
   
zones(12).vertices = [ ...
    xlines(v_PivotRight,v_e,v_PivotLeft-v_PivotRight,v_SplitCenter); ...
    v_PivotRight; ...
    xlines(v_PivotRight,v_RightTop,v_SplitRight,v_LeftTop-v_RightTop); ...
    xlines(v_e,v_RightTop,v_SplitCenter,v_LeftTop-v_RightTop) ];
if abs(dot(v_SplitCenter,v_SplitRight)) == norm(v_SplitCenter)*norm(v_SplitRight)
    zones(12).linearization = @(xypos) [2 ... 
        dot(xypos-segments(2).origin,segments(2).direction)];
else
    xpoint = xlines(v_PivotRight,v_e,v_SplitRight,v_SplitCenter);
    zones(12).linearization = @(xypos) [2 norm( ...
        xlines(xypos,segments(2).origin,xypos-xpoint,segments(2).direction) - ...
        segments(2).origin)];
end

zones(13).vertices = [ ...
    v_PivotLeft; ...
    xlines(v_PivotRight,v_e,v_PivotLeft-v_PivotRight,v_SplitCenter); ...
    xlines(v_e,v_RightTop,v_SplitCenter,v_LeftTop-v_RightTop); ...
    xlines(v_PivotLeft,v_RightTop,v_SplitLeft,v_LeftTop-v_RightTop) ];
if abs(dot(v_SplitCenter,v_SplitLeft)) == norm(v_SplitCenter)*norm(v_SplitLeft)
    zones(13).linearization = @(xypos) [4 ... 
        dot(xypos-segments(4).origin,segments(4).direction)];
else
    xpoint = xlines(v_PivotLeft,v_e,v_SplitLeft,v_SplitCenter);
    zones(13).linearization = @(xypos) [4 norm( ...
        xlines(xypos,segments(4).origin,xypos-xpoint,segments(4).direction) - ...
        segments(4).origin)];
end

%{
% For anyone who is trying to understand the zone definitions: 
% this little loop illustrates the definitions of the zones
for i = 1:length(zones)
    h = patch(zones(i).vertices(:,1),zones(i).vertices(:,2),'r','FaceAlpha',0.5);
    pause;
    if exist('h')
        delete(h);
    end
end
%}

% Now using these zone definitions, we apply inpolygon to determine 
% when the trajectory is within each zone. This zone-occupancy information
% is stored in zoneflags
% Each row of zoneflags corresponds to a position sample in smoothedpos.data,
% whereas each column corresponds to a zone defined in the zone cell array
zoneflags = zeros(size(smoothedpos.data,1),length(zones));
for j = 1:length(zones)
    zoneflags(1:end,j) = logical( inpolygon( ...
        smoothedpos.data(1:end,4),smoothedpos.data(1:end,5), ... 
        zones(j).vertices(:,1),zones(j).vertices(:,2) ) );
end

% Define a size constant for the maximum number of linearizations which
% may possibly apply to the position at any one time. Almost always, 
% this number ranges between 1 and 3, but in case of an amazing 
% coincidence it is possible for this to be as many as 5 (given the
% zones defined above), so just to be safe I will pick this value
MAX_NUM_OVERLAP_ZONES = 5;
linearizations = zeros(size(smoothedpos.data,1),MAX_NUM_OVERLAP_ZONES,2);
linearizations(:,:,2) = -1;

% Compute candidate linearizations over all position samples
for i = 1:size(zoneflags,1)
    currentzones = find(zoneflags(i,:));
    if isempty(currentzones)
        % the position sample does not fall into any of the zones! raise an error!
        error('zone definitions are faulty; this is likely due to bad W-track definition input, or possibly due to a bug');
    end
    for j = 1:numel(currentzones)
        linearizations(i,j,:) = ... 
            zones(currentzones(j)).linearization(smoothedpos.data(i,4:5));
        if linearizations(i,j,2) > segments(linearizations(i,j,1)).length;
            disp(linearizations(i,j,:));
            error('there is a problem with the linearization!');
        end
    end
    % remainder of linearizations(i,:,1) is padded with 0
    % remainder of linearizations(i,:,2) is padded with -1
    % THIS IS VERY IMPORTANT! DO NOT USE NaN values here! Use -1, -23453, or
    % any other negative number. Just not NaN, because NaN always fails
    % equality tests! 
end

% Define an array to contain the consensus linearization
consensus = [ ... 
            zeros(size(smoothedpos.data,1),1) ...
            -ones(size(smoothedpos.data,1),1) ];

% run through the possible linearizations and select from them into the 
% consensus array
while 1
    % retain records of the previous state of the consensus array
    last_consensus = consensus;
    last_linearizations = linearizations;
    % run through position samples one-by-one, forwards and backwards
    idxs = [ 1:size(smoothedpos.data,1) (size(smoothedpos.data,1)-1):(-1):1 ];
    for i = 2:numel(idxs)
        if consensus(idxs(i),1) 
            % if we have assigned a consensus linearization, skip
            continue;
        else
            % STAGE 1: 
            % If the last sample that we checked has a consensus
            % linearization, we use it as a reference for deciding whether
            % to cull invalid linearizations (those which incur impossible
            % jumps in linpos)
            candidates = find(linearizations(idxs(i),:,1));
            if numel(candidates) == 1
                % check: if we have only one linearization candidate, 
                % use it
                consensus(idxs(i),:) = linearizations(idxs(i),1,:);
                continue;
            end
            if consensus(idxs(i-1),1)
                neighbor = consensus(idxs(i-1),:);
                for j = candidates
                    jump = lindist( neighbor, ...
                        shiftdim(linearizations(idxs(i),j,:),1) );
                    if jump > MAXJUMP
                        % wipe out the candidate which violates MAXJUMP
                        linearizations(idxs(i),j,1) = 0;
                        linearizations(idxs(i),j,2) = -1;
                    end
                end
                % and repackage this row of linearizations
                candidates = find(linearizations(idxs(i),:,1));
                if ~isempty(candidates)
                    linearizations(idxs(i),1:numel(candidates),:) = ...
                    linearizations(idxs(i),candidates,:);
                else
                    error(['All candidate linearizations have ' ... 
                    'been disqualified. This may be due to a jump ' ...
                    'in the smoothed position in posparams, or it ' ...
                    'may be due to a yet-uncharacterized bug in ' ...
                    'this function']);
                end
            end
            % STAGE 2: 
            % We check whether the remaining candidate linearizations are
            % within PROXIMITY lindist of each other. If there are pairs
            % which satisfy this proximity criterion, we merge them using
            % linmidpoint.
            while 1
                candidates = find(linearizations(idxs(i),:,1));
                if numel(candidates) == 1
                    consensus(idxs(i),:) = linearizations(idxs(i),1,:);
                    break;
                end
                % Compute pairwise distances between the candidate linearizations.
                % To generate the indices to do this, we construct an array of ones
                % and zeros, where ones are present only for indices corresponding to
                % candidate linearizations, and pairwise comparisons are not 
                % duplicated
                pairs = triu(1 - eye(size(linearizations,2)));
                noncandidates = setdiff(1:size(linearizations,2),candidates);
                pairs(noncandidates,:) = 0;
                pairs(:,noncandidates) = 0;
                [jvec kvec] = find(pairs);
                % compute pairwise distances and find the minimum
                dist = lindist( shiftdim(linearizations(idxs(i),jvec,:),1),...
                    shiftdim(linearizations(idxs(i),kvec,:),1) );
                [mindist, minidx] = min(dist);
                % If the closest pair are within PROXIMITY of each other, and
                % if they lie on the same segment, merge them
                if (mindist < PROXIMITY) && ... 
                ( linearizations(idxs(i),jvec(minidx),1) == ... 
                  linearizations(idxs(i),kvec(minidx),1) )
                    % overwrite one member of this close-proximity pair
                    new = linmidpoint( ...
                        shiftdim(linearizations(idxs(i),jvec(minidx),:),1), ...
                        shiftdim(linearizations(idxs(i),kvec(minidx),:),1) );
                    linearizations(idxs(i),jvec(minidx),:) = new;
                    % and wipe out the other member of the pair
                    linearizations(idxs(i),kvec(minidx),1) = 0;
                    linearizations(idxs(i),kvec(minidx),2) = -1;
                    % and repackage this row of linearizations
                    candidates = find(linearizations(idxs(i),:,1));
                    if ~isempty(candidates)
                        linearizations(idxs(i),1:numel(candidates),:) = ...
                        linearizations(idxs(i),candidates,:);
                    else
                        error(['All candidate linearizations have ' ... 
                        'been disqualified. This may be due to a jump ' ...
                        'in the smoothed position in posparams, or it ' ...
                        'may be due to a yet-uncharacterized bug in ' ...
                        'this function']);
                    end
                else
                    break;
                end
            end
            % STAGE 3: 
            % We do a final check to see whether we have managed to eliminate all
            % but one candidate linearization
            candidates = find(linearizations(idxs(i),:,1));
            if numel(candidates) == 1
                % check: if we have only one linearization candidate, 
                % use it
                consensus(idxs(i),:) = linearizations(idxs(i),1,:);
            end
        end
    end
    % jump out of the while loop if we haven't accomplished anything in this pass
    if all(consensus(:) == last_consensus(:)) && ... 
       all(linearizations(:) == last_linearizations(:))
        break
    end
end

% identify the intervals over which we don't have a consensus linearization
start_idxs = find(diff(consensus(:,1) == 0) > 0) + 1;
if consensus(1,1) == 0
    start_idxs = [1; start_idxs];
end
end_idxs = find(diff(consensus(:,1) == 0) < 0);
if consensus(end,1) == 0
    end_idxs = [end_idxs; size(consensus,1)];
end
unassigned = [start_idxs end_idxs];

% now we let the human operator assign linearization to the position samples
% in these uncertain intervals. to do this, we need to set up some figure
% infrastructure to provide information to the user
main = figure();
axis(posrange);
grid on;
hold on;
trajectory = ... 
    line(smoothedpos.data(:,4),smoothedpos.data(:,5),'Color','c');
line([v_d(1); v_a(1)],[v_d(2); v_a(2)],'Color','k');
line([v_e(1); v_d(1)],[v_e(2); v_d(2)],'Color','k');
line([v_e(1); v_b(1)],[v_e(2); v_b(2)],'Color','k');
line([v_e(1); v_f(1)],[v_e(2); v_f(2)],'Color','k');
line([v_f(1); v_c(1)],[v_f(2); v_c(2)],'Color','k');

global interactivemode;
interactivemode = 1;
global keypress;
keypress = [];
global currentpath;
currentpath = line([NaN],[NaN],'Color','b');
currentpos = line([NaN],[NaN], ...
    'Color','b','MarkerSize',4,'Marker','o','LineWidth',2);
markers = [];
for k = 1:size(linearizations,2)
    markers(k) = line([NaN],[NaN], ...
    'Color','r','MarkerSize',4,'Marker','o','LineWidth',1);
end
currenttime = text(posrange(2)-1,posrange(4)-1,'',... 
    'HorizontalAlignment','right','VerticalAlignment','top');
progress = text(posrange(1)+1,posrange(4)-1,'', ...
    'HorizontalAlignment','left','VerticalAlignment','top');
listcandidates = text(posrange(1)+1,posrange(4)-6,'', ...
    'HorizontalAlignment','left','VerticalAlignment','top');
global animationspeed;
animationspeed = 100;

instructions = [ ...
    'For intervals over which the correct linearization ' ...
    'can not be determined by this function, mark the ' ...
    'correct linearization of each position sample by ' ...
    'hand. The trajectory over each such interval will ' ...
    'shown as an animated loop, which will repeat until ' ...
    'it is interrupted by pressing [spacebar]. Once the ' ...
    'this animated loop is interupted, you can choose ' ...
    'among the candidate linearizations by pressing the ' ...
    'the numbered keys "1", "2", "3" ... (your choice ' ...
    'will change color from red to green. To step to the ' ...
    'next position sample, press "n"; to step to the ' ...
    'previous position sample, press "b". To resume the ' ...
    'animated loop, press [spacebar]. Any linearization ' ...
    'assignments that you make will be stored instantly ' ...
    'in a persistent buffer. When all position samples ' ...
    'within an uncertain interval have been assigned a ' ...
    'linearization, you will be asked whether you want ' ...
    'commit the linearization assignments for the ' ...
    'interval. Once committed, linearization ' ...
    'assignments can not be undone. You can clear all ' ...
    'linearization assignments for an interval by ' ...
    'pressing "c".' ];
popup = msgbox(instructions,'modal');
set(popup,'Resize','on','Units','characters');
popup_contents = get(popup,'Children');
popup_position = get(popup,'Position');
popup_position(3) = size(instructions,2)+5;
popup_position(4) = size(instructions,1)+5;
uiwait(popup);

set(main,'KeyPressFcn',@interceptkeypress);

% now let the human operator patch up each "gap" in the consensus 
% linearization by hand
for i = 1:size(unassigned,1)
    % indices of the position samples within current interval
    idxs = unassigned(i,1):unassigned(i,2);
    % draw the segment of the trajectory corresponding to the 
    % current interval 
    set(currentpath,'XData',smoothedpos.data(idxs,4),...
        'YData',smoothedpos.data(idxs,5));
    % also, draw the trajectory over the preceding and following 
    % 5 seconds (useful for seeing the behavioral context)
    trajidxs = intersect( ... 
        (idxs(1)-ceil(5*FRAMERATE)):(idxs(end)+ceil(5*FRAMERATE)), ...
        1:size(smoothedpos.data,1) );
    set(trajectory,'XData',smoothedpos.data(trajidxs,4),...
        'YData',smoothedpos.data(trajidxs,5));
    % create a vector of integers storing the user's choice of
    % linearization for each position sample in interval.
    assignments = zeros(length(idxs),1);
    % initialize iterator dummy variable for looping through 
    % the current interval 
    j = idxs(1);
    jnext = j+1;
    while 1
        % check to see whether there is a complete sequence of
        % assignments over the entire interval; if so, offer
        % to commit the assignments to consensus
        % update the figure elements to display the
        % current position sample
        set(currentpos,'XData',[smoothedpos.data(j,4)], ...
            'YData',[smoothedpos.data(j,5)]);
        set(currenttime,'String',num2str(smoothedpos.data(j,1)));
        set(progress,'String',[ 'Linearizations ' ... 
            'assigned to ' num2str(nnz(assignments)) ...
            ' out of ' num2str(length(assignments)) ... 
            ' position samples in this interval' ]);
        for k = 1:size(linearizations,2)
            linpos = shiftdim(linearizations(j,k,:),1);
            if linpos(1)
                xypos = lin2cartesian(linpos);
                set(markers(k),'XData',[xypos(1)],'YData',[xypos(2)]);
                if k == assignments(j-idxs(1)+1)
                    set(markers(k),'Color','g');
                else
                    set(markers(k),'Color','r');
                end
            else
                set(markers(k),'XData',[NaN],'YData',[NaN]);
            end
        end
        candidatesstring = {};
        for l = find(linearizations(j,:,1))
            if l == assignments(j-idxs(1)+1)
                candidatesstring{l} = [ '* candidate ' ...
                    num2str(l) ' = segment ' ... 
                    num2str(linearizations(j,l,1)) ... 
                    ' @' num2str(linearizations(j,l,2)) ' *'];
            else
                candidatesstring{l} = [ 'candidate ' ...
                    num2str(l) ' = segment ' ... 
                    num2str(linearizations(j,l,1)) ... 
                    ' @' num2str(linearizations(j,l,2)) ];
            end
        end
        set(listcandidates,'String',candidatesstring);
        if all(assignments)
            set(progress,'String',[ 'Linearizations ' ... 
                'assigned to all ' num2str(length(assignments)) ... 
                ' position samples in this interval; press "w" ' ...
                'to commit' ]);
            if keypress == 'w'
                % perform a lindist test to make sure that none of the
                % selections violate MAXJUMP
                dist = [];
                for l = 1:(length(idxs)-1)
                    jump(l) = lindist( ...
                        shiftdim(linearizations(idxs(l),assignments(l),:),1),...
                        shiftdim(linearizations(idxs(l+1),assignments(l+1),:),1) );
                end
                if any(jump > MAXJUMP)
                    set(progress,'String',[ 'Can not commit linearization ' ...
                        'assignments because some assignments incur ' ...
                        'implausible jumps in linpos; please reassign' ]);
                else
                    break
                end
            end
        end
        if ~interactivemode
            if j == idxs(end)
                jnext = idxs(1);
                interactivemode = 1;
            else
                jnext = j+1;
                if assignments(j-idxs(1)+1) && (abs(jnext-j) == 1) ...
                    && ~assignments(jnext-idxs(1)+1)
                    % find the linearization which is closest to the 
                    % last one
                    candidates = find(linearizations(jnext,:,1));
                    % compute linearized distances of the candidates from the
                    % last assigned position
                    kvec = repmat(assignments(j-idxs(1)+1),numel(candidates),1);
                    dist = [];
                    dist = lindist( ...
                        shiftdim(linearizations(jnext,candidates,:),1),...
                        shiftdim(linearizations(j,kvec,:),1) );
                    [mindist, minidx] = min(dist);
                    % if the nearest candidate is within MAXJUMP, we go ahead
                    % and assign it
                    if (mindist < MAXJUMP)
                        assignments(jnext-idxs(1)+1) = minidx;
                    end
                    if ~assignments(jnext-idxs(1)+1)
                        interactivemode = 1;
                    end
                end
            end
            j = jnext;
            pause(1/animationspeed);
        else
            waitforbuttonpress; 
            if ~isempty(str2num(keypress)) && isscalar(str2num(keypress)) ... 
                   && ismember(str2num(keypress),find(linearizations(j,:,1)))
                selection = str2num(keypress);
                if assignments(j-idxs(1)+1) == selection
                    set(markers(assignments(j-idxs(1)+1)),'Color','r');
                    assignments(j-idxs(1)+1) = 0;
                    set(progress,'String',[ 'Linearizations ' ... 
                        'assigned to ' num2str(nnz(assignments)) ...
                        ' out of ' num2str(length(assignments)) ... 
                        ' position samples in this interval' ]);
                else
                    % all others should be colored red
                    assignments(j-idxs(1)+1) = 0;
                    for k = 1:size(linearizations,2)
                        set(markers(k),'Color','r');
                    end
                    assignments(j-idxs(1)+1) = selection;
                    set(markers(assignments(j-idxs(1)+1)),'Color','g');
                    set(progress,'String',[ 'Linearizations ' ... 
                        'assigned to ' num2str(nnz(assignments)) ...
                        ' out of ' num2str(length(assignments)) ... 
                        ' position samples in this interval' ]);
                end
                % if nothing else, we remain on the same position sample
                % await further instructions from the user
                jnext = j;
            elseif strcmp(keypress,'n')
                if j == idxs(end)
                    jnext = idxs(1);
                else
                    jnext = j+1;
                end
            elseif strcmp(keypress,'b')
                if j == idxs(1)
                    jnext = idxs(end);
                else
                    jnext = j-1;
                end
            elseif strcmp(keypress,'c')
                % clear all assignments for this interval
                assignments(:) = 0;
            end
            if assignments(j-idxs(1)+1) && (abs(jnext-j) == 1) && ...
               ~assignments(jnext-idxs(1)+1)
                % if we don't have an assignment for the next position 
                % sample, we pick one which is closest to the current position
                % sample - assuming that a plausible one can be found
                candidates = find(linearizations(jnext,:,1));
                % compute linearized distances of the candidates from the
                % last assigned position
                kvec = repmat(assignments(j-idxs(1)+1),numel(candidates),1);
                dist = [];
                dist = lindist( ...
                    shiftdim(linearizations(jnext,candidates,:),1),...
                    shiftdim(linearizations(j,kvec,:),1) );
                [mindist, minidx] = min(dist);
                % if the nearest candidate is within MAXJUMP, we go ahead
                % and assign it
                if (mindist < MAXJUMP)
                    assignments(jnext-idxs(1)+1) = minidx;
                end
            end
            j = jnext;
        end
    end
    % once we've finished assigning linearizations to the interval,
    % write the assignments to consensus
    for j = idxs;
        consensus(j,:) = shiftdim( ...
            linearizations(j,assignments(j-idxs(1)+1),:), 1 );
    end
end

delete(main);

% remember that possruct has following fields:
%    time1 raw_x2 raw_y3 smoothed_x4 smoothed_y5 
%    running_speed6 running_direction7 
%    head_direction8 immobility_flag9

linpos.descript = strvcat('linearization based upon', ...
    smoothedpos.descript);
linpos.fields = ... 
    'timestamp1 segment2 linpos3 linspeed4 linheading5 immobility6';
linpos.Wtrack = Wtrack;
linpos.data = NaN(size(smoothedpos.data,1),5);

% inherit timestamps
linpos.data(:,1) = smoothedpos.data(:,1);

% the consensus linearization: segment assignments and lindist along
% each segment
linpos.data(:,2:3) = consensus;

% projection of speed vector onto the segment
speedvecs = zeros(size(smoothedpos.data,1),2);
[speedvecs(:,1) speedvecs(:,2)] = ... 
    pol2cart(smoothedpos.data(:,7),smoothedpos.data(:,6));
for i = 1:size(linpos.data,1)
    linpos.data(i,4) = dot( speedvecs(i,:),... 
    segments(linpos.data(i,2)).direction, 2);
end

% two requirements: 
% movement velocity vector must align with the current track segment
% AND
% head direction vector must align with the current track segment

% compute angle between real head direction and the direction of the current
% track segment. If these angles are within 45degrees of each other, we
% mark the linheading as +1 or -1 (with respect to the direction of the 
% current track segment)
for i = 1:size(linpos.data,1)
    headdir = smoothedpos.data(i,8);
    speeddir = smoothedpos.data(i,7);
    segmentdir = cart2pol( segments(linpos.data(i,2)).direction(1), ...
        segments(linpos.data(i,2)).direction(2) );
    headanglediff = min( ... 
        abs(headdir-segmentdir), 2*pi-abs(headdir-segmentdir) );
    speedanglediff = min( ... 
        abs(speeddir-segmentdir), 2*pi-abs(speeddir-segmentdir) );
    if ( (headanglediff < pi/4) & (speedanglediff < pi/4) )
        linpos.data(i,5) = 1; % head pointing away from segment origin
    elseif ( (headanglediff > 3*pi/4) & (speedanglediff > 3*pi/4) )
        linpos.data(i,5) = -1; % head pointing towards segment origin
    end
end

end % end main function linearizeWtrack


% -----------------------------------------------------------------
function xypos = lin2cartesian(linpos)
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % !!!! VERY IMPORTANT THAT segments BE DECLARED GLOBAL !!!!
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    global segments;
    % This function converts linearized distance along a segment to 
    % position in [x y] coordinates. Note that this function 
    % definition is intimately related to the segments struct array.
    % Given Nx2 array linpos, in which each row is a 2x1 linpos vector,
    % this function returns an Nx2 vector of [x y ] coordinates. Note
    % that this function definition is intimately related to the 
    % segments struct array.
    if (size(linpos,2) ~= 2) 
        error('argument must be Nx2 array');
    elseif any(~ismember(linpos(:,1),1:numel(segments))) 
        error('at least one of the arguments does not represent a valid segment');
    elseif any(linpos(:,2) < 0)
        error('distance along segment must be greater than or equal to zero');
    else 
        for n = 1:numel(segments)
            if any(linpos(find(linpos(:,1) == n),2) > segments(n).length)
                disp(linpos(:,2));
                disp(segments(n).length); 
                error('distance along segment can not exceed the total length of the segment');
            end
        end
    end
    for n = 1:size(linpos,1)  
        xypos(n,:) = segments(linpos(n,1)).origin + ...
        linpos(n,2)*segments(linpos(n,1)).direction;
    end
end % end subfunction lin2cartesian

% -----------------------------------------------------------------
function d = lindist(a,b)
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % !!!! VERY IMPORTANT THAT segments BE DECLARED GLOBAL !!!!
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    global segments;
    % Given Nx2 arrays a and b, in which each row is a 2x1 linpos vector,
    % this function returns an Nx1 vector of linear distances between
    % corresponding rows of a and b (distances constrained along track
    % segments, sort of like a "city-block" distance metric). This is
    % very specific to the segment definitions, so any changes made
    % in how the segments are defined must be reflected here
    if (size(a,2) ~= 2) || (size(b,2) ~= 2)
        error('argument must be Nx2 array');
    elseif any(~ismember(a(:,1),1:numel(segments))) || ... 
        any(~ismember(b(:,1),1:numel(segments)))
        error('at least one of the arguments does not represent a valid segment');
    elseif any(a(:,2) < 0) || any(b(:,2) < 0)
        error('distance along segment must be greater than or equal to zero');
    else 
        for n = 1:numel(segments)
            if any( a(find(a(:,1) == n),2) > segments(n).length ) || ...
                any( b(find(b(:,1) == n),2) > segments(n).length )
                error('distance along segment can not exceed the total length of the segment');
            end
        end
    end
    for n = 1:size(a,1)  
        switch a(n,1)
        % easy case is when both linpos arguments share segment
        case b(n,1)
            d(n) = abs(a(n,2) - b(n,2));
        % otherwise, jump into a big LUT - this could be optimized, but
        % at expense of code clarity. anyway, on most calls to this 
        % function, flow will not go beyond the initial if block
        case 1
            switch b(n,1)
            case 2
                d(n) = a(n,2) + segments(2).length - b(n,2);
            case 3
                d(n) = a(n,2) + segments(2).length + b(n,2);
            case 4
                d(n) = a(n,2) + segments(2).length + b(n,2);
            case 5
                d(n) = a(n,2) + segments(2).length + segments(4).length + b(n,2);
            end
        case 2
            switch b(n,1)
            case 1
                d(n) = segments(2).length - a(n,2) + b(n,2);
            case 3
                d(n) = a(n,2) + b(n,2);
            case 4
                d(n) = a(n,2) + b(n,2);
            case 5
                d(n) = a(n,2) + segments(4).length + b(n,2);
            end
        case 3
            switch b(n,1)
            case 2
                d(n) = a(n,2) + b(n,2);
            case 4
                d(n) = a(n,2) + b(n,2);
            case 1
                d(n) = a(n,2) + segments(2).length + b(n,2);
            case 5
                d(n) = a(n,2) + segments(4).length + b(n,2);
            end
        case 4
            switch b(n,1)
            case 5
                d(n) = segments(4).length - a(n,2) + b(n,2);
            case 3
                d(n) = a(n,2) + b(n,2);
            case 2
                d(n) = a(n,2) + b(n,2);
            case 1
                d(n) = a(n,2) + segments(2).length + b(n,2);
            end
        case 5
            switch b(n,1)
            case 4
                d(n) = a(n,2) + segments(4).length - b(n,2);
            case 3
                d(n) = a(n,2) + segments(4).length + b(n,2);
            case 2
                d(n) = a(n,2) + segments(4).length + b(n,2);
            case 1
                d(n) = a(n,2) + segments(4).length + segments(2).length + b(n,2);
            end
        end 
    end
end % end subfunction lindist


% -----------------------------------------------------------------
function m = linmidpoint(a,b)
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % !!!! VERY IMPORTANT THAT segments BE DECLARED GLOBAL !!!!
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    global segments;
    % Given Nx2 arrays a and b, in which each row is a 2x1 linpos vector,
    % this function returns an Nx2 vector of "midpoints" between
    % corresponding rows of a and b (constrained along track segments, 
    % in "city-block" metric space). This is very specific to the segment 
    % definitions, so any changes made in how the segments are defined 
    % must be reflected here
    if (size(a,2) ~= 2) || (size(b,2) ~= 2)
        error('argument must be Nx2 array');
    elseif any(~ismember(a(:,1),1:numel(segments))) || ... 
        any(~ismember(b(:,1),1:numel(segments)))
        error('at least one of the arguments does not represent a valid segment');
    elseif any(a(:,2) < 0) || any(b(:,2) < 0)
        error('distance along segment must be greater than or equal to zero');
    else 
        for n = 1:numel(segments)
            if any( a(find(a(:,1) == n),2) > segments(n).length ) || ...
                any( b(find(b(:,1) == n),2) > segments(n).length )
                error('distance along segment can not exceed the total length of the segment');
            end
        end
    end
    for n = 1:size(a,1)    
        % easy case is when both linpos arguments share segment
        if a(n,1) == b(n,1)
            m(n,:) = [ a(n,1) (a(n,2)+b(n,2))/2 ];
        else
            % first, compute one-half the linearized distance between
            % a(n,:) and b(n,:)
            md = 0.5*lindist(a(n,:),b(n,:));
            % now jump into a big LUT which specifies how many vertices
            % intervene between a(n,:) and b(n,:) and determines which 
            % segment to place the midpoint
            switch a(n,1)
            case 1
                switch b(n,1)
                case 2
                    if md > a(n,2)
                        m(n,:) = [2, segments(2).length - (md - a(n,2)) ];
                    else
                        m(n,:) = [1, a(n,2) - md ];
                    end 
                case 3
                    if md > a(n,2) + segments(2).length
                        m(n,:) = [3, md - segments(2).length - a(n,2) ];
                    elseif md > a(n,2)
                        m(n,:) = [2, segments(2).length - (md - a(n,2)) ];
                    else
                        m(n,:) = [1, a(n,2) - md ];
                    end
                case 4
                    if md > a(n,2) + segments(2).length
                        m(n,:) = [4, md - segments(2).length - a(n,2) ];
                    elseif md > a(n,2)
                        m(n,:) = [2, segments(2).length - (md - a(n,2)) ];
                    else
                        m(n,:) = [1, a(n,2) - md ];
                    end
                case 5
                    if md > a(n,2) + segments(2).length + segments(4).length
                        m(n,:) = [5, md - segments(4).length - segments(2).length - a(n,2) ];
                    elseif md > a(n,2) + segments(2).length 
                        m(n,:) = [4, md - segments(2).length - a(n,2) ];
                    elseif md > a(n,2)
                        m(n,:) = [2, segments(2).length - (md - a(n,2)) ];
                    else
                        m(n,:) = [1, a(n,2) - md ];
                    end
                end
            case 2
                switch b(n,1)
                case 1
                    if md > segments(2).length - a(n,2)
                        m(n,:) = [1, md - (segments(2).length - a(n,2)) ];
                    else
                        m(n,:) = [2, a(n,2) + md ];
                    end
                case 3
                    if md > a(n,2)
                        m(n,:) = [3, md - a(n,2) ];
                    else
                        m(n,:) = [2, a(n,2) - md ];
                    end
                case 4
                    if md > a(n,2)
                        m(n,:) = [4, md - a(n,2) ];
                    else
                        m(n,:) = [2, a(n,2) - md ];
                    end
                case 5
                    if md > a(n,2) + segments(4).length
                        m(n,:) = [5, md - segments(4).length - a(n,2) ];
                    elseif md > a(n,2)
                        m(n,:) = [4, md - a(n,2) ];
                    else
                        m(n,:) = [2, a(n,2) - md ];
                    end
                end
            case 3
                switch b(n,1)
                case 2
                    if md > a(n,2)
                        m(n,:) = [2, md - a(n,2) ];
                    else
                        m(n,:) = [3, a(n,2) - md ];
                    end
                case 4
                    if md > a(n,2)
                        m(n,:) = [4, md - a(n,2) ];
                    else 
                        m(n,:) = [3, a(n,2) - md ];
                    end
                case 1
                    if md > a(n,2) + segments(2).length
                        m(n,:) = [1, md - segments(2).length - a(n,2) ];
                    elseif md > a(n,2)
                        m(n,:) = [2, md - a(n,2) ];
                    else
                        m(n,:) = [3, a(n,2) - md ];
                    end
                case 5
                    if md > a(n,2) + segments(4).length
                        m(n,:) = [5, md - segments(4).length - a(n,2) ];
                    elseif md > a(n,2)
                        m(n,:) = [4, md - a(n,2) ];
                    else
                        m(n,:) = [3, a(n,2) - md ];
                    end
                end
            case 4
                switch b(n,1)
                case 5
                    if md > segments(4).length - a(n,2)
                        m(n,:) = [5, md - (segments(4).length - a(n,2)) ];
                    else
                        m(n,:) = [4, a(n,2) + md ];
                    end
                case 3
                    if md > a(n,2)
                        m(n,:) = [3, md - a(n,2) ];
                    else
                        m(n,:) = [4, a(n,2) - md ];
                    end
                case 2
                    if md > a(n,2)
                        m(n,:) = [2, md - a(n,2) ];
                    else
                        m(n,:) = [4, a(n,2) - md ];
                    end
                case 1
                    if md > a(n,2) + segments(2).length
                        m(n,:) = [1, md - segments(2).length - a(n,2) ];
                    elseif md > a(n,2)
                        m(n,:) = [2, md - a(n,2) ];
                    else
                        m(n,:) = [4, a(n,2) - md ];
                    end
                end
            case 5
                switch b(n,1)
                case 4
                    if md > a(n,2)
                        m(n,:) = [4, segments(4).length - (md - a(n,2)) ];
                    else
                        m(n,:) = [5, a(n,2) - md ];
                    end 
                case 3
                    if md > a(n,2) + segments(4).length
                        m(n,:) = [3, md - segments(4).length - a(n,2) ];
                    elseif md > a(n,2)
                        m(n,:) = [4, segments(4).length - (md - a(n,2)) ];
                    else
                        m(n,:) = [5, a(n,2) - md ];
                    end
                case 2
                    if md > a(n,2) + segments(4).length
                        m(n,:) = [2, md - segments(4).length - a(n,2) ];
                    elseif md > a(n,2)
                        m(n,:) = [4, segments(4).length - (md - a(n,2)) ];
                    else
                        m(n,:) = [5, a(n,2) - md ];
                    end
                case 1
                    if md > a(n,2) + segments(4).length + segments(2).length
                        m(n,:) = [1, md - segments(2).length - segments(4).length - a(n,2) ];
                    elseif md > a(n,2) + segments(4).length
                        m(n,:) = [2, md - segments(4).length - a(n,2) ];
                    elseif md > a(n,2)
                        m(n,:) = [4, segments(4).length - (md - a(n,2)) ];
                    else
                        m(n,:) = [5, a(n,2) - md ];
                    end
                end
            end
        end
        if m(n,2) > segments(m(n,1)).length % final check
            error('there is a bug in linmidpoint');
        end
    end
end % end subfunction linmidpoint


% -----------------------------------------------------------------
function u = perp(v)
    % for a 1x2 vector v, returns a vector u which is perpendicular
    if size(v) ~= [1 2]
        error('perp requires a 1x2 row vector as an argument')
    end
    u = [ -v(2) v(1) ];
end % end subfunction perp


% -----------------------------------------------------------------
function p3 = xlines(p1,p2,v_a,v_b)
    % Line-Line Intersection
    % p3 = xlines(p1,p2,v_a,v_b)
    % define a subfunction for finding the intersection of two lines
    % L1 and L2 defined in 2-space, where L1 has direction v_a and
    % contains point p1, and L2 has direction v_b and contains point
    % p2. All arguments are row vectors of the form [x y].
    if ( size(p1) ~= [1 2] | ...
         size(p2) ~= [1 2] | ...
         size(v_a) ~= [1 2] | ...
         size(v_b) ~= [1 2])
        error('xlines requires 1x2 row vectors as arguments');
    end
    A = [v_a; -v_b];
    if abs(det(A*A')) > eps % check that the lines are not parallel
        parallelflag = 0;
        d = A'\(p2 - p1)'; % least squares pseudo-inverse
        % we measure from p1 and from p2 and average the two values
        p3 = 0.5*(p1 + d(1)*v_a) + 0.5*(p2 + d(2)*v_b);
    else
        p3 = [];
    end 
end % end subfunction xlines

% -----------------------------------------------------------------
function interceptkeypress(obj,eventdata)
    global interactivemode;
    global keypress;
    global animationspeed;
    % this callback processes user keyboard input
    keypress = get(obj,'CurrentCharacter');
    if isempty(keypress)
        keypress = '';
    elseif ~interactivemode && strcmp(keypress,' ')
        interactivemode = 1;
    elseif interactivemode && strcmp(keypress,' ');
        interactivemode = 0;
    elseif strcmp(keypress,'+') && (animationspeed < 1000)
        animationspeed = animationspeed*1.5;
    elseif strcmp(keypress,'-') && (animationspeed > 0.5)
        animationspeed = animationspeed/1.5;
    end
end % end callback function interceptkeypress
