function sss_presort_pca(spikes, mode, nch)
% Shantanu: Apr 2011

if (~ismember(mode, {'xy','xyz'})), error('Unknown mode.');  end;


if nch==1,
    % Use 1st 4 PCs
    for i=1:4
        cmd=sprintf('data(:,i)=spikes.pcadata.pc%d;',i); eval(cmd);
    end
end

if nch==2,
    % Use 1st 2 PCs of the two channels
    data(:,1)=spikes.pcadata.pc1(:,1);
    data(:,2)=spikes.pcadata.pc1(:,2);
    data(:,3)=spikes.pcadata.pc2(:,1);
    data(:,4)=spikes.pcadata.pc2(:,2);
end

if nch==3
     % Use 1st PC of 3 channels + an additional 2nd PC
    for i=1:3
        data(:,i)=spikes.pcadata.pc1(:,i);
    end
    data(:,4)=spikes.pcadata.pc2(:,1);
end


if nch>=4,
    % Use 1st PC of 4 channels
    for i=1:4
        data(:,i)=spikes.pcadata.pc1(:,i);
    end
end

% Make the plot.
figure;  hfig = gcf;  hax = gca;    hold on;
redimscreen70s;
%orient(gcf,'landscape');
set(gcf, 'PaperPositionMode', 'auto');
set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');

inds{1}=1:size(data,1);
if (strcmp(mode, 'xyz'))
    hndl = plot3(data(:,1), data(:,2), data(:,3),'.','MarkerSize',2);
else
    hndl = plot(data(:,1), data(:,2), '.','MarkerSize',2);
    set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
end

h_pts = hndl;
ssghandle = hndl;

if (length(find(h_pts)) == 1), set(h_pts, 'MarkerSize', 2, 'Marker', '.');  end;
if (~strcmp(mode,'xyz')),  uistack(h_pts, 'top');  end;

hold off;

% Make the figure's spike-sorting gui (SSG) object.
ssg.mode = mode;
ssg.ss_object = spikes;
ssg.ss_assigns = 1*(size(data,1));
ssg.group_indices = inds; %1*(size(data,1));
ssg.group_handles = ssghandle;
ssg.xchoice = 'PC';     ssg.xparam1 = '1';      ssg.xcontrol = [];
ssg.ychoice = 'PC';     ssg.yparam1 = '2';      ssg.ycontrol = [];
ssg.zchoice = 'PC';     ssg.zparam1 = '3';      ssg.zcontrol = [];
guidata(gcf, ssg);

% Colorize each axis for easy identification.
set(hax, 'XColor', [0.8 0 0.5], 'YColor', [0.1 0.7 0.3], 'ZColor', [0 0.4 0.7]);

% Make the axis labels & callback-ify them.
textprops = {'FontSize', 14, 'FontWeight', 'bold', 'Rotation', 0};
hx = xlabel([ssg.xchoice ssg.xparam1]);  set(hx, 'ButtonDownFcn', {@make_control, 'x'}, textprops{:});
hy = ylabel([ssg.ychoice ssg.yparam1]);  set(hy, 'ButtonDownFcn', {@make_control, 'y'}, textprops{:});
if (strcmp(mode,'xyz'))
    hz = zlabel([ssg.zchoice ssg.zparam1]);  set(hz, 'ButtonDownFcn', {@make_control, 'z'}, textprops{:});
    cameratoolbar('ResetCamera');  
    cameratoolbar('SetMode', 'orbit');
	set(gcf, 'Renderer', 'OpenGL');
    grid on;
else
    cameratoolbar('SetMode', 'nomode');
	set(gcf, 'Renderer', 'zbuffer');
end
axis equal;

set(hfig, 'DeleteFcn', @delete_function);
set(hfig, 'ButtonDownFcn', {@make_density, hax});
figure(hfig);  % bring figure to top


function delete_function(hObject, event)
% Need to delete any associated control GUIs.
ssg = guidata(hObject);
delete([ssg.xcontrol, ssg.ycontrol, ssg.zcontrol]);


function make_control(hObject, event, controlaxis)
% Raises the associated axis control if it exists and creates it if it does not.
% Note that the axis control is responsible for registering its handle with the
% the guidata of this axis upon its creation/deletion.
ssg = guidata(hObject);
controlhandle = ssg.([controlaxis 'control']);
if (~isempty(controlhandle))
    figure(controlhandle);
else
    ssg_featureselect(gca, controlaxis);
end

function make_density(myfig, event, myaxes)
% Replots the current data as a density using histxy or hist3d.
ssg = guidata(myfig);
if(strcmp(get(myfig, 'SelectionType'), 'open'))
    [az,el] = view;   T = view;
    properties2d = {'XLim', 'YLim', 'DataAspectRatio'};
	properties3d = cat(2, properties2d, {'Zlim', 'CameraViewAngle'});
    properties2d = cat(1, properties2d, get(gca,properties2d));
	properties3d = cat(1, properties3d, get(gca,properties3d));
	
    xdata = [];    ydata = [];   zdata = [];
    for clust = 1:length(ssg.group_handles)
        xdata = [xdata get(ssg.group_handles(clust), 'XData')];
        ydata = [ydata get(ssg.group_handles(clust), 'YData')];
        if (strcmp(ssg.mode,'xyz'))
            zdata = [zdata get(ssg.group_handles(clust), 'ZData')];
        end
    end
    figure; hdens = gca;
    if (strcmp(ssg.mode,'xyz'))
        view(T);  set(gca,properties3d{:});
		levels = logspace(-3,0,10);
		hist3d([xdata', ydata', zdata'], [], [], 20, [], 0);  
		set(hdens, 'Box', 'on');  grid on;
    else
		histxy(xdata,ydata,250,1);
        set(gca,properties2d{:},'Color',[0 0 0.508]);
        HGlogmenu(findobj(gca,'Type','Image'));
    end
end

