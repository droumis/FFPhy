function registeredpos = registerpos_Utrack(rawpos,template)
%
%   registerpos_Utrack(rawpos,Utrack_template);
% 
% This requires the MATLAB Image Processing toolbox

registeredpos = rawpos;
registeredpos.descript = 'position data registered to template';
registeredpos = rmfield(registeredpos,{'sampleframe','source'});
registeredpos.data(:,1) = rawpos.data(:,1);

main = figure();
image(template.controlpoints_guide);
set(gca,'DataAspectRatio',[1 1 1],'YDir','reverse','Visible','off');
popup = msgbox('click on control points in this sequence; press [space] to confirm selection','modal');
popup_contents = get(popup,'Children');
uiwait(popup);

close(gcf);
% fit axes to the size of the image frame
image([0.5 size(rawpos.sampleframe,2)-0.5], ...
    [0.5 size(rawpos.sampleframe,1)-0.5], ...
    rawpos.sampleframe,'AlphaData',1); 
set(gca,'DataAspectRatio',[1 1 1],'YDir','normal','Box','on');

% define a UserData struct which packages relevant info
UserData.controlpoints = ginput(size(template.controlpoints,1));

% draw the segments of the track and collect object handles
% in the array linehdl
UserData.linehdl = [];
UserData.linehdl(1) = line( ...
    [UserData.controlpoints(1,1); UserData.controlpoints(2,1)], ... 
    [UserData.controlpoints(1,2); UserData.controlpoints(2,2)], ...
    'Color','r','Marker','o','MarkerSize',4, ...
    'MarkerFaceColor','r','HitTest','off');
UserData.linehdl(2) = line( ...
    [UserData.controlpoints(2,1); UserData.controlpoints(6,1)], ... 
    [UserData.controlpoints(2,2); UserData.controlpoints(6,2)], ...
    'Color','r','Marker','o','MarkerSize',4, ...
    'MarkerFaceColor','r','HitTest','off');
UserData.linehdl(3) = line( ...
    [UserData.controlpoints(6,1); UserData.controlpoints(7,1)], ... 
    [UserData.controlpoints(6,2); UserData.controlpoints(7,2)], ...
    'Color','r','Marker','o','MarkerSize',4, ...
    'MarkerFaceColor','r','HitTest','off');
UserData.linehdl(4) = line( ...
    [UserData.controlpoints(7,1); UserData.controlpoints(3,1)], ... 
    [UserData.controlpoints(7,2); UserData.controlpoints(3,2)], ...
    'Color','r','Marker','o','MarkerSize',4, ...
    'MarkerFaceColor','r','HitTest','off');
UserData.linehdl(5) = line( ...
    [UserData.controlpoints(3,1); UserData.controlpoints(4,1)], ... 
    [UserData.controlpoints(3,2); UserData.controlpoints(4,2)], ...
    'Color','r','Marker','o','MarkerSize',4, ...
    'MarkerFaceColor','r','HitTest','off');
UserData.linehdl(6) = line( ...
    [UserData.controlpoints(4,1); UserData.controlpoints(8,1)], ... 
    [UserData.controlpoints(4,2); UserData.controlpoints(8,2)], ...
    'Color','r','Marker','o','MarkerSize',4, ...
    'MarkerFaceColor','r','HitTest','off');
UserData.linehdl(7) = line( ...
    [UserData.controlpoints(8,1); UserData.controlpoints(5,1)], ... 
    [UserData.controlpoints(8,2); UserData.controlpoints(5,2)], ...
    'Color','r','Marker','o','MarkerSize',4, ...
    'MarkerFaceColor','r','HitTest','off');
UserData.linehdl(8) = line( ...
    [UserData.controlpoints(5,1); UserData.controlpoints(1,1)], ... 
    [UserData.controlpoints(5,2); UserData.controlpoints(1,2)], ...
    'Color','r','Marker','o','MarkerSize',4, ...
    'MarkerFaceColor','r','HitTest','off');
% last location of the mouse cursor
UserData.lastxy = [];
% currently-selected vertex
UserData.selectedvertex = 0;
% flag for whether the user is finished with input
UserData.finished = 0;
% GUI callback hooks
set(main,'KeyPressFcn',@interceptkeypress, ...
'WindowButtonDownFcn',@mousebuttondown,'UserData',UserData);

% now let the user adjust the controlpoints as desired until
% the [spacebar] is pressed to commmit the controlpoints
while 1
    waitforbuttonpress;
    UserData = get(main,'UserData');
    if UserData.finished
        break;
    end
end
close(gcf);

tform = cp2tform(UserData.controlpoints,template.controlpoints,'projective');
registeredpos.data(:,2:3) = tform.forward_fcn(rawpos.data(:,2:3),tform);
try
    registeredpos.data(:,4:5) = tform.forward_fcn(rawpos.data(:,4:5),tform);
catch
    %pass
end

registeredpos.registeredframe = imtransform(rawpos.sampleframe,tform);
registeredpos.template = template;
registeredpos.transform = tform;

end

%----------------------------------------------------------------------
function mousebuttondown(obj,eventdata)
    % grab a local copy of UserData
    UserData = get(obj,'UserData');
    % define a proximity within which a mouse click
    % captures the nearest vertex
    axeshandle = get(obj,'CurrentAxes');
    set(obj,'Units','pixels');
    set(axeshandle,'Units','pixels');
    axes_xlim = get(axeshandle,'XLim'); % axis coords (cm)
    axes_position = get(axeshandle,'Position'); % axes in screen pixels
    axes_width = axes_position(3); % pixels
    scale_x = (axes_xlim(2) - axes_xlim(1))/axes_width; % cm/pixel
    HITRADIUS = 10*scale_x; % what distance does 10 pixels correspond to?
    % reassign callbacks
    set(obj,'WindowButtonMotionFcn',@mousedrag);
    set(obj,'WindowButtonDownFcn','');
    set(obj,'WindowButtonUpFcn',@mousebuttonup);
    % figure out where the current cursor is
    UserData.lastxy = getxycursorpos(obj);
    dist = hypot( ...
        UserData.controlpoints(:,1) - UserData.lastxy(1), ...
        UserData.controlpoints(:,2) - UserData.lastxy(2) );
    [mindist, minidx] = min(dist);
    if mindist < HITRADIUS
        UserData.selectedvertex = minidx;
    else
        UserData.selectedvertex = 0;
    end
    % copy UserData back to main
    set(obj,'UserData',UserData);
end % end callback mousebuttondown


%----------------------------------------------------------------------
function mousebuttonup(obj,eventdata)
    % grab a local copy of UserData
    UserData = get(obj,'UserData');
    set(obj,'WindowButtonDownFcn',@mousebuttondown);
    set(obj,'WindowButtonMotionFcn','');
    UserData.lastxy = [];
    % copy UserData back to main
    set(obj,'UserData',UserData);
end % end callback mousebuttonup

%----------------------------------------------------------------------
function mousedrag(obj,eventdata)
    % grab a local copy of UserData
    UserData = get(obj,'UserData');
    axes_xlim = get(get(obj,'CurrentAxes'),'XLim');
    axes_ylim = get(get(obj,'CurrentAxes'),'YLim');
    set(obj,'WindowButtonDownFcn','');
    set(obj,'WindowButtonUpFcn',@mousebuttonup);
    % adjust controlpoints by the amount that CurrentPoint changed
    newxy = getxycursorpos(obj);
    if ~isempty(UserData.lastxy) && UserData.selectedvertex && ...
    (newxy(1) > axes_xlim(1)) && (newxy(1) < axes_xlim(2)) && ...
    (newxy(2) > axes_ylim(1)) && (newxy(2) < axes_ylim(2))
        UserData.controlpoints(UserData.selectedvertex,:) = ...
        UserData.controlpoints(UserData.selectedvertex,:) + (newxy - UserData.lastxy);
        set(UserData.linehdl(1), ...
            'XData',[UserData.controlpoints(1,1); UserData.controlpoints(2,1)], ... 
            'YData',[UserData.controlpoints(1,2); UserData.controlpoints(2,2)]);
        set(UserData.linehdl(2), ...
            'XData',[UserData.controlpoints(2,1); UserData.controlpoints(6,1)], ... 
            'YData',[UserData.controlpoints(2,2); UserData.controlpoints(6,2)]);
        set(UserData.linehdl(3), ...
            'XData',[UserData.controlpoints(6,1); UserData.controlpoints(7,1)], ... 
            'YData',[UserData.controlpoints(6,2); UserData.controlpoints(7,2)]);
        set(UserData.linehdl(4), ...
            'XData',[UserData.controlpoints(7,1); UserData.controlpoints(3,1)], ... 
            'YData',[UserData.controlpoints(7,2); UserData.controlpoints(3,2)]);
        set(UserData.linehdl(5), ...
            'XData',[UserData.controlpoints(3,1); UserData.controlpoints(4,1)], ... 
            'YData',[UserData.controlpoints(3,2); UserData.controlpoints(4,2)]);
        set(UserData.linehdl(6), ...
            'XData',[UserData.controlpoints(4,1); UserData.controlpoints(8,1)], ... 
            'YData',[UserData.controlpoints(4,2); UserData.controlpoints(8,2)]);
        set(UserData.linehdl(7), ...
            'XData',[UserData.controlpoints(8,1); UserData.controlpoints(5,1)], ... 
            'YData',[UserData.controlpoints(8,2); UserData.controlpoints(5,2)]);
        set(UserData.linehdl(8), ...
            'XData',[UserData.controlpoints(5,1); UserData.controlpoints(1,1)], ... 
            'YData',[UserData.controlpoints(5,2); UserData.controlpoints(1,2)]);
        UserData.lastxy = newxy;
    end
    % copy UserData back to main
    set(obj,'UserData',UserData);
end % end callback mousedrag

%----------------------------------------------------------------------
function xy = getxycursorpos(obj)
% query the figure to convert cursor position in the figure 
% window to [x y] coordinates in the coordinate frame of the
% plot axes
    cp = get(get(obj,'CurrentAxes'),'CurrentPoint');
    xy(1) = cp(1,1);
    xy(2) = cp(1,2);
end % end subfunction getxycursorpos

%----------------------------------------------------------------------
function interceptkeypress(obj,eventdata)
    % grab a local copy of UserData
    UserData = get(obj,'UserData');
    keypress = get(obj,'CurrentCharacter');
    if strcmp(keypress,' ')
        UserData.finished = 1; 
    end
    % copy UserData back to main
    set(obj,'UserData',UserData);
end % end callback interceptkeypress
