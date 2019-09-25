function Wtrack = defineWtrack(rawpos)
%
%   Wtrack = defineWtrack(rawpos)
%
% rawpos is a rawpos structure for a run session. It is 
% important that there be no gaps in the rawpos, and that 
% all timestamps are in monotonically increasing order.
%
% The output, Wtrack, is a struct with two fields:
%
%   posrange is [xmin xmax ymin ymax], the 
%   rectangle which bounds all position values
%
%   Wtrack.vertices is a 5x2 array which contains 
%   the coordinates of the vertices of the Wtrack. Each row
%   contains the [x y] coordinates of a vertex. Vertices 
%   are listed in the order shown below:
%
%          6-----5-----4     
%          |     |     |     
%          |     |     |     
%          |     |     |     
%          |     |     |     
%          |     |     |     
%          3     2     1     
%
%

% Error checking on the rawpos input
if any(isnan(rawpos.data(:,4:5)))
    error('rawpos data contains NaN values; you need to first interpolate over these');
elseif any(diff(rawpos.data(:,1)) < 0)
    error('rawpos data contains out-of-order timestamps');
end

% the output struct
global Wtrack;
Wtrack.descript = 'Wtrack definition';
Wtrack.subject = rawpos.subject;
Wtrack.day = rawpos.day;
Wtrack.session = rawpos.session;
Wtrack.environment = rawpos.environment;

% Get the user to specify the vertices of the W track
posrange = [ ...
    min(rawpos.data(:,4))-20, ...
    max(rawpos.data(:,4))+20, ...
    min(rawpos.data(:,5))-20, ...
    max(rawpos.data(:,5))+20 ];
instructions = [ ...
    'Using the graphical input crosshairs,              ' ; ...
    'click on the vertices of the W track               ' ; ...
    'in the sequence indicated  below:                  ' ; ...
    '                                                   ' ; ...
    '          6-----5-----4                            ' ; ...
    '          |     |     |                            ' ; ...
    '          |     |     |                            ' ; ...
    '          |     |     |                            ' ; ...
    '          |     |     |                            ' ; ...
    '          3     2     1                            ' ; ...
    '                                                   ' ; ...
    '(Note that the reference track above               ' ; ...
    'may not match the orientation of your              ' ; ...
    'position data. Be sure that your                   ' ; ...
    'vertex selections respect the                      ' ; ...
    'chirality of the track.)                           ' ; ...
    '                                                   ' ; ...
    'You can drag the vertices to adjust                ' ; ...
    'them. When you are satisfied, press                ' ; ...
    '[spacebar] to commit.                              ' ];
popup = msgbox(instructions,'modal');
popup_contents = get(popup,'Children');
set(get(popup_contents(1),'Children'),'FontName','fixedwidth');
set(popup,'Resize','on','Units','characters');
popup_position = get(popup,'Position');
popup_position(3) = size(instructions,2)+5;
popup_position(4) = size(instructions,1)+5;
set(popup,'Position',popup_position);
uiwait(popup);

% draw a main figure window showing trajectory
main = figure();
line(rawpos.data(:,4),rawpos.data(:,5),'Color','r');
axis(posrange);
grid on;
set(main,'KeyPressFcn',@interceptkeypress,... 
    'WindowButtonDownFcn',@mousebuttondown);

% we construct a state struct which stores various 
% state variables in its fields
% let the user pick out the vertices
state.vertices = ginput(6);
% draw the segments of the track and collect object handles
% in the array linehdl
state.linehdl = [];
state.linehdl(1) = line( ...
    [state.vertices(4,1); state.vertices(1,1)], ... 
    [state.vertices(4,2); state.vertices(1,2)], ...
    'Color','k','Marker','s','MarkerSize',4, ...
    'MarkerFaceColor','k','HitTest','off');
state.linehdl(2) = line( ...
    [state.vertices(4,1); state.vertices(5,1)], ...
    [state.vertices(4,2); state.vertices(5,2)], ...
    'Color','k','Marker','s','MarkerSize',4, ...
    'MarkerFaceColor','k','HitTest','off');
state.linehdl(3) = line( ... 
    [state.vertices(5,1); state.vertices(2,1)], ...
    [state.vertices(5,2); state.vertices(2,2)], ...
    'Color','k','Marker','s','MarkerSize',4, ...
    'MarkerFaceColor','k','HitTest','off');
state.linehdl(4) = line( ...
    [state.vertices(6,1); state.vertices(5,1)], ...
    [state.vertices(6,2); state.vertices(5,2)], ...
    'Color','k','Marker','s','MarkerSize',4, ...
    'MarkerFaceColor','k','HitTest','off');
state.linehdl(5) = line( ...
    [state.vertices(6,1); state.vertices(3,1)], ...
    [state.vertices(6,2); state.vertices(3,2)], ...
    'Color','k','Marker','s','MarkerSize',4, ...
    'MarkerFaceColor','k','HitTest','off');
% last location of the mouse cursor
state.lastxy = [];
% currently-selected vertex
state.selectedvertex = 0;
% flag for whether the user is finished with input
state.finished = 0;
% save this state data in the UserData property of main
set(main,'UserData',state);

message = text(posrange(2)-1,posrange(4)-1, ...
    '','HorizontalAlignment','right','VerticalAlignment','top');

% now let the user adjust the vertices as desired until
% the [spacebar] is pressed to commmit the vertices
while 1
    waitforbuttonpress;
    set(message,'String','');
    state = get(main,'UserData');
    if state.finished
        MAXDIFF = 2; % maximum difference in the lengths of arms (cm)
        if abs( norm(state.vertices(6,:) - state.vertices(5,:)) - ... 
                norm(state.vertices(4,:) - state.vertices(5,:)) ) ...
               > MAXDIFF
            set(message,'String', ...
                'Track arms are of uneven length. Please fix.');
            state.finished = 0;
            set(main,'UserData',state);
        else
            break;
        end
    end
end

% write the final vertex definitions to the output
Wtrack.vertices = state.vertices;

% Construct a large collection of vectors for partitioning the trajectory
% into zones with different linearization schemes
v_a = Wtrack.vertices(1,:);
v_b = Wtrack.vertices(2,:);
v_c = Wtrack.vertices(3,:);
v_d = Wtrack.vertices(4,:);
v_e = Wtrack.vertices(5,:);
v_f = Wtrack.vertices(6,:);

% define vectors in the directions along the track segments; note that
% this also provides information on the lengths of the track segments
v_1 = v_a - v_d;
v_2 = v_d - v_e;
v_3 = v_b - v_e;
v_4 = v_f - v_e;
v_5 = v_c - v_f;

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
% compute segment definitions.
Wtrack.segments = [];

Wtrack.segments(1).origin = v_d;
Wtrack.segments(1).direction = v_1/norm(v_1);
Wtrack.segments(1).length = Inf;

Wtrack.segments(2).origin = v_e;
Wtrack.segments(2).direction = v_2/norm(v_2);
Wtrack.segments(2).length = norm(v_2);

Wtrack.segments(3).origin = v_e;
Wtrack.segments(3).direction = v_3/norm(v_3);
Wtrack.segments(3).length = Inf;

Wtrack.segments(4).origin = v_e;
Wtrack.segments(4).direction = v_4/norm(v_4);
Wtrack.segments(4).length = norm(v_4);

Wtrack.segments(5).origin = v_f;
Wtrack.segments(5).direction = v_5/norm(v_5);
Wtrack.segments(5).length = Inf;

delete(main);
end % end main function
    
%----------------------------------------------------------------------
function mousebuttondown(obj,eventdata)
    % define a proximity within which a mouse click
    % captures the nearest vertex
    axeshandle = get(obj,'CurrentAxes');
    set(obj,'Units','pixels');
    set(axeshandle,'Units','pixels');
    axes_xlim = get(axeshandle,'XLim'); % axis coords (cm)
    axes_position = get(axeshandle,'Position');
    axes_width = axes_position(3); % pixels
    scale_x = (axes_xlim(2) - axes_xlim(1))/axes_width; % cm/pixel
    HITRADIUS = 3*scale_x; % what does 3 pixels correspond to?
    % grab a local copy of the state struct
    state = get(obj,'UserData');
    set(obj,'WindowButtonMotionFcn',@mousedrag);
    set(obj,'WindowButtonDownFcn','');
    set(obj,'WindowButtonUpFcn',@mousebuttonup);
    % figure out where the current cursor is
    state.lastxy = convertcursorpos(obj);
    dist = hypot( ...
        state.vertices(:,1) - state.lastxy(1), ...
        state.vertices(:,2) - state.lastxy(2) );
    [mindist, minidx] = min(dist);
    if mindist < HITRADIUS
        state.selectedvertex = minidx;
    else
        state.selectedvertex = 0;
    end
    % copy the local copy of state back to the parent object
    set(obj,'UserData',state);
end % end callback mousebuttondown


%----------------------------------------------------------------------
function mousebuttonup(obj,eventdata)
    % grab a local copy of the state struct
    state = get(obj,'UserData');
    set(obj,'WindowButtonDownFcn',@mousebuttondown);
    set(obj,'WindowButtonMotionFcn','');
    state.lastxy = [];
    % copy the local copy of state back to the parent object
    set(obj,'UserData',state);
end % end callback mousebuttonup

%----------------------------------------------------------------------
function mousedrag(obj,eventdata)
    axes_xlim = get(get(obj,'CurrentAxes'),'XLim');
    axes_ylim = get(get(obj,'CurrentAxes'),'YLim');
    % grab a local copy of the state struct
    state = get(obj,'UserData');
    set(obj,'WindowButtonDownFcn','');
    set(obj,'WindowButtonUpFcn',@mousebuttonup);
    % adjust vertices by the amount that CurrentPoint changed
    newxy = convertcursorpos(obj);
    if ~isempty(state.lastxy) && state.selectedvertex && ...
       (newxy(1) > axes_xlim(1)) && ...
       (newxy(1) < axes_xlim(2)) && ...
       (newxy(2) > axes_ylim(1)) && ...
       (newxy(2) < axes_ylim(2))
        state.vertices(state.selectedvertex,:) = ... 
            state.vertices(state.selectedvertex,:) + ...
            (newxy - state.lastxy);
        set(state.linehdl(1), ...
            'XData',[state.vertices(4,1); state.vertices(1,1)], ... 
            'YData',[state.vertices(4,2); state.vertices(1,2)]);
        set(state.linehdl(2), ...
            'XData',[state.vertices(4,1); state.vertices(5,1)], ...
            'YData',[state.vertices(4,2); state.vertices(5,2)]);
        set(state.linehdl(3), ...
            'XData',[state.vertices(5,1); state.vertices(2,1)], ...
            'YData',[state.vertices(5,2); state.vertices(2,2)]);
        set(state.linehdl(4), ...
            'XData',[state.vertices(6,1); state.vertices(5,1)], ...
            'YData',[state.vertices(6,2); state.vertices(5,2)]);
        set(state.linehdl(5), ...
            'XData',[state.vertices(6,1); state.vertices(3,1)], ...
            'YData',[state.vertices(6,2); state.vertices(3,2)]);
        state.lastxy = newxy;
    end
    % copy the local copy of state back to the parent object
    set(obj,'UserData',state);
end % end callback mousedrag

%----------------------------------------------------------------------
function xy = convertcursorpos(obj)
% query the figure to convert cursor position in the figure 
% window to [x y] coordinates in the coordinate frame of the
% plot axes
    cursorpos = get(obj,'CurrentPoint');
    axeshandle = get(obj,'CurrentAxes');
    set(obj,'Units','pixels');
    set(axeshandle,'Units','pixels');
    axes_xlim = get(axeshandle,'XLim'); % cm
    axes_ylim = get(axeshandle,'YLim'); % cm
    axes_position = get(axeshandle,'Position'); % pixels
    axes_left = axes_position(1);
    axes_bottom = axes_position(2);
    axes_width = axes_position(3);
    axes_height = axes_position(4);
    scale_x = (axes_xlim(2) - axes_xlim(1))/axes_width; % cm/pixel
    scale_y = (axes_ylim(2) - axes_ylim(1))/axes_height; % cm/pixel
    xy(1) = axes_xlim(1) + scale_x*(cursorpos(1) - axes_left);
    xy(2) = axes_ylim(1) + scale_y*(cursorpos(2) - axes_bottom);
end % end subfunction convertcursorpos

%----------------------------------------------------------------------
function interceptkeypress(obj,eventdata)
    % grab a local copy of the state struct
    state = get(obj,'UserData');
    keypress = get(obj,'CurrentCharacter');
    if strcmp(keypress,' ')
        state.finished = 1; 
    end
    % copy the local copy of state back to the parent object
    set(obj,'UserData',state);
end % end callback interceptkeypress
