function sss_presort_amp(spikes, mode, nch)
% 2d or 3d Amplitude Plots
%  Shantanu: Apr2011

if (~ismember(mode, {'xy','xyz'})), error('Unknown mode.');  end;


% if nch==1
%     disp('Single Channel Only - Plotting Amplitude against time')
%     ampl = abs(max(spikes.waveforms,[],2)) + abs(min(spikes.waveforms,[],2));
%     spktimes = spikes.fstimes./1000; % in sec
%     
%     afig=figure; hold on;
%     set(gcf, 'PaperPositionMode', 'auto');
%     set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
%     ylabel('Amplitude (uV)' ,'FontSize', 14, 'FontWeight', 'bold');
%     xlabel('Time (sec)','FontSize', 14, 'FontWeight', 'bold');
%     
%     hndl = plot(spktimes, ampl, '.');
%     title('Amplitude vs Time')
%     hold off;
%     
%     return
% end

channel_lth = size(spikes.waveforms,2)/nch;
nSpikes = size(spikes.waveforms,1);

% Calculate Amplitude for mulit-channel data
if ~isfield(spikes,'ampl')
    
    for ch=1:nch
        if isfield(spikes,'waveforms_ch1')
            cmd=sprintf('ampl(:,ch)=abs(max(spikes.waveforms_ch%d,[],2)) + abs(min(spikes.waveforms_ch%d,[],2));',ch,ch); eval(cmd);
        else
            w = spikes.waveforms(:,channel_lth*(ch-1)+1:channel_lth*ch);
            ampl(:,ch) = abs(max(w,[],2)) + abs(min(w,[],2));
        end
    end    
else
    ampl=spikes.ampl;
end


% Make the plot.
figure;  hfig = gcf;  hax = gca;    hold on;
redimscreen70s;
%orient(gcf,'landscape');
set(gcf, 'PaperPositionMode', 'auto');
set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');


if nch==1,
    disp('Single Channel: 2nd dimension (X axis) is Time ');
    ampl(:,2)=ampl(:,1);
    ampl(:,1)=spikes.fstimes./1000;
end

if nch==2,
    disp('Stereotrode: 3rd dimension (Z axis) is Time ');
    ampl(:,3)=spikes.fstimes./1000;
end

data=ampl;
spikes.ampl=ampl;
inds{1}=1:size(data,1);
if nch==1,
    disp('Cannot do 3-d plot for 1-channel data. Plotting 2-d instead')
    mode='xy';
end
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
ssg.group_indices = inds;
ssg.group_handles = ssghandle;
ssg.xchoice = 'Amp';     ssg.xparam1 = '1';      ssg.xcontrol = [];
ssg.ychoice = 'Amp';     ssg.yparam1 = '2';      ssg.ycontrol = [];
ssg.zchoice = 'Amp';     ssg.zparam1 = '3';      ssg.zcontrol = [];
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

% % Set up the legend (unless there are too many clusters, in which case just show noise).
% clustnames = cellstr(num2str(show));
% if (exist('h_out', 'var')),  h_pts = [h_out, h_pts];   clustnames = cat(1, {'Outliers'}, clustnames);  end;
%
% %if (numclusts < 12  &&  numclusts > 1)
% if (length(show) < 12  &&  length(show) > 1)    
%     hleg = legend(h_pts, clustnames, 0);
% elseif (exist('h_out', 'var')),
%     hleg = legend(h_pts(1), clustnames{1}, 0);
% end
% if (strcmp(mode,'xyz') && exist('hleg'))
%     legpos = get(hleg, 'Position');
%     legpos(2) = 0.78;
%     set(hleg, 'Position', legpos);
% end

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
    %     ssg_ampfeatureselect(gca, controlaxis);
    ssg_ampfeatureselect(gca, controlaxis);
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

