function zonesstruct = definezones(smoothedpos)
%
%   zones = definezones(smoothedpos);
%

% Error checking on the posparams input
if any(isnan(smoothedpos.data(:,4:5)))
    error('posparams data contains NaN values; you need to first interpolate over these');
elseif any(diff(smoothedpos.data(:,1)) < 0)
    error('posparams data contains out-of-order timestamps');
end

instructions = 'define zones of interest';
popup = helpdlg(instructions,'');
uiwait(popup);

global main;
main = figure();
set(main,'KeyPressFcn',@interceptkeypress,... 
    'WindowButtonDownFcn',@mousebuttondown);

% two tools for the user to define a zone of interest
%
%   a circle tool: draw a near-circular polygon centered on the
%   selected point. when the user clicks on the workspace, a dialog box 
%   will ask for a radius; if the default value of 0 is accepted, then 
%   the user will need to use the mouse to adjust the radius manually
%
%   a polygon tool: draw a polygon. clicks on the workspace will be
%   interpreted as vertices of the polygon. a middle-button click will
%   close the polygon
%
% when a zone definition is completed, ask user for a description string
% create radio button group for selecting zones
global circle_tool;
global polygon_tool;
tool_panel = uibuttongroup('Position',[0 0.95 1 0.05], ...
    'Visible','off','HandleVisibility','callback',...
    'SelectionChangeFcn',@select_tool);
circle_tool = uicontrol('Style','Radio','String','circle',...
    'Parent',tool_panel,'HandleVisibility','off', ...
    'Units','normalized','Position',[0 0 0.5 1]);
polygon_tool = uicontrol('Style','Radio','String','manual polygon',...
    'Parent',tool_panel,'HandleVisibility','off', ...
    'Units','normalized','Position',[0.5 0 0.5 1]);
set(tool_panel,'SelectedObject',circle_tool);  
set(tool_panel,'Visible','on');

% draw position trajectory in the main window
axes('XLim',[min(smoothedpos.data(:,4))-20 max(smoothedpos.data(:,4))+20], ...
    'YLim',[min(smoothedpos.data(:,5))-20 max(smoothedpos.data(:,5))+20], ...
    'XGrid','on','YGrid','on','ActivePositionProperty','OuterPosition', ...
    'PlotBoxAspectRatioMode','auto','DataAspectRatio',[1 1 1], ...
    'Units','normalized');
line(smoothedpos.data(:,4),smoothedpos.data(:,5),'Color','r');

% we construct a state struct which stores various state variables
% the zone definitions. note that numel(state.zones) is the number of defined zones
state.zones = struct([]);
% last location of the mouse cursor
state.lastxy = [];
% which tool is currently selected in the radiobutton panel?
state.currenttool = circle_tool;
% flag for whether the user is finished (triggers quit)
state.finished = 0;
% also construct a vector of handles to the zone patch objects
state.patches = [];
% save this state data in the UserData property of the main gui object
set(main,'UserData',state);

% now let the user adjust the vertices as desired until
% the [spacebar] is pressed to commmit the vertices
while 1
    waitforbuttonpress;
    state = get(main,'UserData');
    if state.finished
        break;
    end
end

% write the final zone definitions to the output
zonesstruct = state.zones;

delete(main);
end % end main function

%----------------------------------------------------------------------
function select_tool(obj,eventdata)
    % grab a local copy of the state struct
    global main;
    state = get(main,'UserData');
    % update current tool
    state.currenttool = eventdata.NewValue;
    % copy the local copy of state back to the parent object
    set(main,'UserData',state);
end % end select_tool callback


%----------------------------------------------------------------------
function mousebuttondown(obj,eventdata)
    global circle_tool;
    global polygon_tool;

    % grab a local copy of the state struct
    state = get(obj,'UserData');

    % if we don't have any defined zones, or if the last zone definition has 
    % been finalized, then propose a new zone definition
    if ~numel(state.zones) || ~isempty(state.zones(end).descript)
        if state.currenttool == circle_tool
            % if currenttool is circle, create a new zone with a single vertex at
            % the current mouse location (center of the circle)
            state.lastxy = convertcursorpos(obj);
            state.zones(end+1).vertices = state.lastxy;
            state.zones(end).descript = [];
            state.patches(end+1) = patch( ...
                'Faces',1:size(state.zones(end).vertices,1), ...
                'Vertices',state.zones(end).vertices, ...
                'FaceColor', [0 0.9 0],'FaceAlpha',0.2,'EdgeColor',[0 1 0]);
            tmp = inputdlg({'enter radius'},'define radius of circular zone',1,{'0'});
            radius = str2num(tmp{1});
            if radius
                % if user specifies a non-zero radius, we use that radius
                % to lock the zone definition, first ask the user to complete the zone definition 
                % by providing a descript string
                tmp = inputdlg({'enter zone name'},'define descript string for zone',1,{' '});
                tmp = tmp{1};
                if isempty(tmp)
                    state.zones(end).descript = ' '; % non-empty string!
                else
                    state.zones(end).descript = tmp;
                end
                center = state.lastxy;
                state.zones(end).vertices = makecircle(center,radius);
                set(state.patches(end), ...
                    'Faces',1:size(state.zones(end).vertices,1), ...
                    'Vertices',state.zones(end).vertices);
                % we're done with this zone definition!
            end
        elseif state.currenttool == polygon_tool
            % if currenttool is polygon, create a new zone with duplicate vertex
            % at the current mouse location
            state.lastxy = convertcursorpos(obj);
            state.zones(end+1).vertices = [state.lastxy; state.lastxy];
            state.zones(end).descript = [];
            state.patches(end+1) = patch( ...
                'Faces',1:size(state.zones(end).vertices,1), ...
                'Vertices',state.zones(end).vertices, ...
                'FaceColor', [0 0.9 0],'FaceAlpha',0.2,'EdgeColor',[0 1 0]);
        else
            error('currenttool does not have a valid value');
        end
    else % ...or update the current (incomplete) zone definition
        if state.currenttool == circle_tool
            % ask the user to complete the zone definition by providing a descript string
            tmp = inputdlg({'enter zone name'},'define descript string for zone',1,{' '});
            tmp = tmp{1};
            if isempty(tmp)
                state.zones(end).descript = ' '; % non-empty string!
            else
                state.zones(end).descript = tmp;
            end
            % set the radius of the circle
            center = state.lastxy;
            radius = norm(center - convertcursorpos(obj));
            state.zones(end).vertices = makecircle(center,radius);
            set(state.patches(end), ...
                'Faces',1:size(state.zones(end).vertices,1), ...
                'Vertices',state.zones(end).vertices);
            % we're done with this zone definition!
        elseif state.currenttool == polygon_tool
            % did the user click the middle button? this signals that the polygon is complete
            if strcmp(get(obj,'SelectionType'),'extend')
                % ask the user to complete the zone definition by providing a descript string
                tmp = inputdlg({'enter zone name'},'define descript string for zone',1,{' '});
                tmp = tmp{1};
                if isempty(tmp)
                    state.zones(end).descript = ' '; % non-empty string!
                else
                    state.zones(end).descript = tmp;
                end
                % we're done with this zone definition!
            end
            state.lastxy = state.zones(end).vertices(end,:);
            % if currenttool is polygon, add vertex to the current zone definition
            state.zones(end).vertices = [ state.zones(end).vertices; state.lastxy ];
            set(state.patches(end), ...
                'Faces',1:size(state.zones(end).vertices,1), ...
                'Vertices',state.zones(end).vertices);
        else
            error('currenttool does not have a valid value');
        end
    end
    set(obj,'CurrentObject',state.patches(end));
    set(obj,'WindowButtonDownFcn','');
    set(obj,'WindowButtonUpFcn',@mousebuttonup);
    set(obj,'WindowButtonMotionFcn',@mousedrag);
    % copy the local copy of state back to the parent object
    set(obj,'UserData',state);
end % end callback mousebuttondown


%----------------------------------------------------------------------
function mousebuttonup(obj,eventdata)
    % grab a local copy of the state struct
    state = get(obj,'UserData');
    set(obj,'WindowButtonUpFcn','');
    set(obj,'WindowButtonDownFcn',@mousebuttondown);
    set(obj,'WindowButtonMotionFcn',@mousedrag);
    % copy the local copy of state back to the parent object
    set(obj,'UserData',state);
end % end callback mousebuttonup

%----------------------------------------------------------------------
function mousedrag(obj,eventdata)
    global circle_tool;
    global polygon_tool;

    % grab a local copy of the state struct
    state = get(obj,'UserData');
    set(obj,'WindowButtonDownFcn',@mousebuttondown);
    set(obj,'WindowButtonUpFcn',@mousebuttonup);
  
    % we only apply this action if a zone exists, and that zone has not been
    % finalized yet
    if numel(state.zones) && isempty(state.zones(end).descript)
        axes_xlim = get(get(obj,'CurrentAxes'),'XLim');
        axes_ylim = get(get(obj,'CurrentAxes'),'YLim');
        % adjust vertices by the amount that CurrentPoint changed
        newxy = convertcursorpos(obj);
        if (newxy(1) > axes_xlim(1)) && ...
            (newxy(1) < axes_xlim(2)) && ...
            (newxy(2) > axes_ylim(1)) && ...
            (newxy(2) < axes_ylim(2))
            % update vertices of the current zone under definition
            if state.currenttool == circle_tool
                center = state.lastxy;
                radius = norm(center - newxy);
                state.zones(end).vertices = makecircle(center,radius);
                set(state.patches(end), ...
                    'Faces',1:size(state.zones(end).vertices,1), ...
                    'Vertices',state.zones(end).vertices);
            elseif state.currenttool == polygon_tool
                state.zones(end).vertices(end,:) = newxy;
                set(state.patches(end), ...
                    'Faces',1:size(state.zones(end).vertices,1), ...
                    'Vertices',state.zones(end).vertices);
            else
                error('currenttool does not have a valid value');
            end
        end
    end
    % copy the local copy of state back to the parent object
    set(obj,'UserData',state);
end % end callback mousedrag

%----------------------------------------------------------------------
function circlevertices = makecircle(center,radius)
% generate a set of vertices for a polygon which approximates a circle
% radius here has units of the plotted position values (i.e. cm)
    numpoints = ceil(radius*4);
    if numpoints < 8
        numpoints = 8;
    end
    t = linspace(0,2*pi,numpoints)';
    r = ones(numpoints,1)*radius;
    [x, y] = pol2cart(t,r);
    circlevertices = [ x+center(1) y+center(2) ];
end % end subfunction circlevertices

%----------------------------------------------------------------------
function xy = convertcursorpos(obj)
% query the figure to convert cursor position in the figure 
% window to [x y] coordinates in the coordinate frame of the
% plot axes
    cursorpos = get(get(obj,'CurrentAxes'),'CurrentPoint');
    xy(1) = cursorpos(1,1);
    xy(2) = cursorpos(1,2);
end % end subfunction convertcursorpos

%----------------------------------------------------------------------
function interceptkeypress(obj,eventdata)

    % grab a local copy of the state struct
    state = get(obj,'UserData');
    keypress = get(obj,'CurrentCharacter');

    switch keypress

    % if user presses [spacebar], check whether all zone definitions are complete
    case ' '
        % if the last zone definition is incomplete, then it will not have a descript field
        % [] as its descript field
        if numel(state.zones) && isempty(state.zones(end).descript)
            msgbox('complete the current zone definition first','modal');
        else
            state.finished = 1;
        end

    % delete the last zone definition
    case 'd'
        if numel(state.zones) > 1
            state.zones = state.zones(1:end-1);
            delete(state.patches(end));
            state.patches(end) = [];
        elseif numel(state.zones) == 1
            state.zones = struct([]);
            delete(state.patches(end));
            state.patches(end) = [];
        end

    end % end switch block

    % copy the local copy of state back to the parent object
    set(obj,'UserData',state);
end % end callback interceptkeypress


