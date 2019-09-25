function cb_logmenu(handle, event)
%CB_LOGMENU        Callback for HGlogmenu.

%%%%%%%%%%%%%%%%%%%%%%%% Determine toggle state %%%%%%%%%%%%%%%%%%%%%%
userdata=get(handle, 'UserData');
switch(get(handle,'Checked')),
    case 'on',  toggle = 'off';  scale = 'linear';
    case 'off', toggle = 'on';   scale = 'log';
end
set(handle, 'Checked',  toggle); 

%%%%%%%%%%%%%%%%%%%%%%%%%% Change the display %%%%%%%%%%%%%%%%%%%%%%%%
switch(userdata.limits(1)),
    case 'Y',
        set(userdata.axes, 'YScale', scale);
    case 'Z',
        set(userdata.axes, 'ZScale', scale);
    case 'C',  % this one is more work
		% first swap limits
		temp = get(userdata.axes, userdata.limits);  
		set(userdata.axes, userdata.limits, userdata.backlims);
		userdata.backlims = temp;   % remember last limits setting ...
		% then swap data
        temp = get(userdata.imageobj, 'CData');
        set(userdata.imageobj, 'CData', userdata.backdata);
        userdata.backdata = temp;
		% then refresh colorbar if one exists
		h = find_colorbar(userdata.axes);
		if (~isempty(h)),
			cbar_pos = get(h,'Position');     axes_pos = get(userdata.axes,'Position');
			if (cbar_pos(1) > sum(axes_pos([1,3]))), loc = 'vert'; else, loc = 'horiz'; end;
			colorbar(loc, 'peer', userdata.axes); 
		end;
		
    otherwise, error('HGlogmenu error: invalid mode.');
end

set(handle, 'UserData', userdata); 