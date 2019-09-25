function h = HGlogmenu(h)
%HGlogmenu         Adds a context menu for toggling logarithmic scaling.
%   HGlogmenu creates a context menu that allows the user to toggle
%   logarithmic scaling for image, surface or line data, depending on the
%   contents of the current axes.
%  
%   If the current axes contain an image, the context menu is associated
%   with that object (or the topmost one if multiple images are present).
%   Pixel intensities toggle between their values when HGlogmenu is called
%   and the log (base 10) of those values (raw data values <= 0 result in
%   NaNs).  The log data is precomputed when the menu is first added to an
%   image to speed up subsequent switching.
%
%   If no image is present but the axes contain surface objects, the
%   context menu is associated with the axes.  The menu then toggles
%   linear vs logarithmic scaling on the z-axis.
%
%   If no image or surface objects are present but the axes contain line
%   objects, the context menu is associated with the axes.  In this case,
%   the menu toggles logarithmic scaling on the y-axis.
%
%   HGlogmenu(HANDLE) associates the menu with the object specified by
%   HANDLE. 
%
%   H = HGlogmenu(...) returns a handle to the object (image or axes)
%   associated with the new context menu.

state = warning;  warning off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin > 1)
    error('Invalid syntax.');
elseif (nargin == 1)
    switch(get(h,'Type')),
        case 'image',    limits = 'CLim';
        case 'surface',  limits = 'ZLim';
        case 'line',     limits = 'YLim';
        otherwise, 
            error(sprintf('HGlogmenu is undefined for objects of type %s.', get(h,'Type')));
    end
else
    childs = get(gca, 'Children');
    img = findobj(childs, 'Type', 'image');
    if (~isempty(img)),
        h = img(1);  limits = 'CLim';
    elseif (~isempty(findobj(childs, 'Type', 'surface'))),
        h = gca;     limits = 'ZLim';
    elseif (~isempty(findobj(childs, 'Type', 'line'))),
        h = gca;     limits = 'YLim';
    else
        error('The current axes do not contain the appropriate objects.');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% Define the Callback %%%%%%%%%%%%%%%%%%%%%%%%
cmenu = get(h, 'UIContextMenu');
if (isempty(cmenu)), cmenu = uicontextmenu; set(h, 'UIContextMenu', cmenu);  end;
logmenu = findobj(cmenu, 'Tag', 'logmenu');
if (isempty(logmenu)),   % if menu doesn't exist yet ...
	logmenu = uimenu(cmenu, 'Tag', 'logmenu', 'Checked', 'off', ...
                            'Label', 'Log scaling', 'Callback', @CB_logmenu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Store UserData %%%%%%%%%%%%%%%%%%%%%%%%%%
switch(limits)
    case {'YLim','ZLim'}
        userdata.axes     = h;
    case 'CLim',
        userdata.axes     = get(h, 'Parent');
        userdata.imageobj = h;
        userdata.backdata = log10(get(h, 'CData'));
            userdata.backdata(get(h,'CData') <= 0) = NaN;   % blank out imaginary/-INF values
        userdata.backlims = [min(userdata.backdata(:)), max(userdata.backdata(:))];  % default tight color scaling on log plot    
end
userdata.limits = limits;
set(logmenu, 'UserData', userdata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cleanup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargout == 0), clear h;  end

warning(state);
