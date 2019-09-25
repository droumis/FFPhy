function linpos = linearizepos_Utrack(smoothedpos,Utrack_template)
%
%   linpos = linearizepos_Utrack(smoothedpos,Utracktemplate)
%
% This function performs an automated linearization given smoothed 
% position data which is registered to the coordinate axes in a 
% particular way. THERE ARE HARD-CODED CONSTANTS IN THIS FUNCTION!!!!
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
if any(isnan(smoothedpos.data(:,2:3)))
    error('posparams data contains NaN values; you need to first interpolate over these');
elseif any(diff(smoothedpos.data(:,1)) <= 0)
    error('posparams data contains out-of-order or duplicate timestamps');
elseif any(hypot( ...
       diff(smoothedpos.data(:,2)),diff(smoothedpos.data(:,3))) > MAXJUMP)
    error('posparams data must not contain any large jumps in position');
end

% Linpos coordinates convention:
%
% segment 1: lindist is measured from vertex "c"
% segment 2: lindist is measured from vertex "d"
% segment 3: lindist is measured from vertex "d"
%
%          d--2--c     
%          |     |     
%          |     |     
%          3     1     
%          |     |     
%          |     |     
%          b     a     
%
segments = Utrack_template.segments;

zones(1).vertices = [ -22 -22; 22 -22; 22 125; -22 125 ];
zones(1).linearization = @(xypos) ...
    [ 1, 147-xypos(2) ];

zones(2).vertices = [ -22 125; 22 125; 22 169; ];
zones(2).linearization = @(xypos) ...
    [ 1, norm( xlines([0 147],[-22 125],[0 -1],xypos-[-22 125]) - [0 147] ) ];

zones(3).vertices = [ -22 125; 22 169; -66 169 ];
zones(3).linearization = @(xypos) ...
    [ 2, norm( xlines([-44 147],[-22 125],[1 0],xypos-[-22 125]) - [-44 147] ) ];
    
zones(4).vertices = [-66 125; -22 125; -66 169 ];
zones(4).linearization = @(xypos) ...
    [ 3, norm( xlines([-44 147],[-22 125],[0 -1],xypos-[-22 125]) - [-44 147] ) ];

zones(5).vertices = [ -66 -22; -22 -22; -22 125; -66 125 ];
zones(5).linearization = @(xypos) ...
    [ 3, 147-xypos(2) ];

%{
% For anyone who is trying to understand the zone definitions: 
% this little loop illustrates the definitions of the zones
h = figure;
line(smoothedpos.data(:,2),smoothedpos.data(:,3),'Color','k');
set(gca,'DataAspectRatio',[1 1 1]);
for i = 1:length(zones)
    patch(zones(i).vertices(:,1),zones(i).vertices(:,2),'r','FaceAlpha',0.5);
    pause(1);
end
delete(h);
%}

% remember that smoothedpos has following fields:
%    time1 smoothed_x2 smoothed_y3 ... 
linpos.descript = 'linearized position data';
linpos.subject = smoothedpos.subject;
linpos.environment = smoothedpos.environment;
linpos.day = smoothedpos.day;
linpos.session = smoothedpos.session;
linpos.tstart = smoothedpos.tstart;
linpos.tend = smoothedpos.tend;
linpos.fields = ... 
    'timestamp[1] segment[2] linpos[3]';
linpos.template = Utrack_template;
linpos.data = NaN(size(smoothedpos.data,1),3);

% inherit timestamps
linpos.data(:,1) = smoothedpos.data(:,1);

% Now using zone definitions, we apply inpolygon to determine 
% when the trajectory is within each zone.
% Each row of zoneflags corresponds to a position sample in smoothedpos.data,
% whereas each column corresponds to a zone defined in the zone cell array
for j = 1:length(zones)
    idx = find(inpolygon( ...
        smoothedpos.data(:,2),smoothedpos.data(:,3), ... 
        zones(j).vertices(:,1),zones(j).vertices(:,2)));
    for i = 1:numel(idx)
        xypos = smoothedpos.data(idx(i),2:3);
        linpos.data(idx(i),2:3) = zones(j).linearization(xypos);
    end
end

end % end main function linearizeWtrack

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

