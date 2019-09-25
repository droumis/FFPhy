function ssgtest(spikes, assignments, show, mode)
% temporary script to translate the SSG_DATABROWSE functions into a GUI.

%   Last Modified By: sbm on Wed Aug 17 18:29:18 2005

if (~ismember(mode, {'xy','xyz'})), error('Unknown mode.');  end;

if (isempty(assignments))
    if (isfield(spikes, 'hierarchy') && isfield(spikes.hierarchy, 'assigns'))
        assignments = spikes.hierarchy.assigns;
    elseif (isfield(spikes, 'overcluster'))
        assignments = spikes.overcluster.assigns;
    else
        assignments = ones(size(spikes.spiketimes));
    end
end
clusters = unique(assignments); numclusts = length(clusters); 

if (isempty(show)),  show = clusters;  end;
% if (isfield(spikes, 'overcluster') && all(ismember(assignments, spikes.overcluster.assigns)))
%     cmap = spikes.overcluster.colors;
% else

    
    cmap=jet(length(clusters));
    %cmap = jet(2*length(show));

%end

%%numclusts = length(clusters);
%numclusts = length(show);   %%%% MY CHANGE: SHANTANU FEB 2006

show(show == 0) = [];


%%%%% PRE-CALCULATED AND STORED PCAS

data=spikes.pca.scores; %%% TAKE 1ST PC from sj_pcasvd

%data=spikes.pca.scores;
% SVD the data because we want to use the PCs as default axes.
% [pca.scores,pca.u,pca.s,pca.v] = pcasvd(spikes.waveforms);
% spikes.pca = pca;
% data = pca.scores;

% Make the plot.
figure;  hfig = gcf;  hax = gca;    hold on;
redimscreen70s;
%orient(gcf,'landscape');
set(gcf, 'PaperPositionMode', 'auto');
set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
    
for clu = 1:numclusts
    sel = find(assignments == clusters(clu));
	inds{clu} = sel;
    if (strcmp(mode, 'xyz'))
        hndl = plot3(data(sel,1), data(sel,2), data(sel,3),'.');
    else
        hndl = plot(data(sel,1), data(sel,2), '.');
        set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    end
	
	if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.1);
	elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:),'MarkerSize',2);
    else                                     set(hndl, 'Color', Clgy, 'MarkerSize',0.1);
    %else                                     set(hndl, 'Color', 'w', 'MarkerSize',0.1);    
    end;
	
	if (clusters(clu) == 0),                 h_out = hndl;
	elseif (ismember(clusters(clu), show)),  h_pts(clusters(clu)) = hndl;
	end
	
	ssghandle(clu) = hndl;
end
h_pts = h_pts(show);

if (length(find(h_pts)) == 1), set(h_pts, 'MarkerSize', 2, 'Marker', '.');  end;
if (~strcmp(mode,'xyz')),  uistack(h_pts, 'top');  end;

hold off;

% Make the figure's spike-sorting gui (SSG) object.
ssg.mode = mode;
ssg.ss_object = spikes;
ssg.ss_assigns = assignments;
ssg.group_indices = inds;   ssg.group_handles = ssghandle;
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

% Set up the legend (unless there are too many clusters, in which case just show noise).
clustnames = cellstr(num2str(show));
if (exist('h_out', 'var')),  h_pts = [h_out, h_pts];   clustnames = cat(1, {'Outliers'}, clustnames);  end;

%if (numclusts < 12  &&  numclusts > 1)
if (length(show) < 12  &&  length(show) > 1)    %% SHANTANU: FEB 2006: CHANGE
    hleg = legend(h_pts, clustnames, 0);    
elseif (exist('h_out', 'var')),
    hleg = legend(h_pts(1), clustnames{1}, 0);
end
if (strcmp(mode,'xyz') && exist('hleg'))
    legpos = get(hleg, 'Position');
    legpos(2) = 0.78;
    set(hleg, 'Position', legpos);
end

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

