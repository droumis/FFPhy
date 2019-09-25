function sss_twoparamplot(spikes, Param1, Param2, assignments, show, nch)
% temporary script to translate the SSG_DATABROWSE functions into a GUI.

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
if (nargin < 4), nch=4; end

ParamNames = {'Time','Amp1','Amp2','Amp3','Amp4','HtCh','WidCh','X-Pos','Y-Pos'};

cmap=jet(length(clusters));
%cmap = hot(2*length(show))

show(show == 0) = [];
%%% TO NOT SHOW OUTLIERS
%clusters(clusters == 0) = []; 
numclusts = length(clusters); 

if isfield(spikes,'matclustparams'),
    data = spikes.matclustparams; datanames=[];
    for i=1:9, datanames{i}=spikes.matclustparamnames{i}; end
    %disp('Using these parameters');
    %disp(spikes.matclustparamnames');
else
    disp('ERROR! No Matlcust Parameters - Using Amplitude')
    for i=1:nch     
        cmd=sprintf('ampl(:,i)=abs(max(spikes.waveforms_ch%d,[],2)) + abs(min(spikes.waveforms_ch%d,[],2));',i,i); eval(cmd);
    end
    data=ampl;
end

% Make the plot.
figure;  hfig = gcf;  hax = gca;    hold on;
redimscreen70s;
%orient(gcf,'landscape');
set(gcf, 'PaperPositionMode', 'auto');
set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
    
for clu = 1:numclusts
    sel = find(assignments == clusters(clu));
	inds{clu} = sel;
    hndl = plot(data(sel,Param1), data(sel,Param2), '.','MarkerSize',2);
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    	
	if (clusters(clu) == 0),                 set(hndl, 'Color', Clgy, 'Marker', '.','MarkerSize',0.1);
	elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:),'MarkerSize',6);
    else                                     set(hndl, 'Color', [0 0 0],'MarkerSize',0.1);
    end;
	
	if (clusters(clu) == 0),                 h_out = hndl;
	elseif (ismember(clusters(clu), show)),  h_pts(clusters(clu)) = hndl;
	end
	ssghandle(clu) = hndl;
end

xlabel(ParamNames{Param1});
ylabel(ParamNames{Param2});

h_pts = h_pts(show);
if (length(find(h_pts)) == 1), set(h_pts, 'MarkerSize', 2, 'Marker', '.');  end;

hold off;
%axis equal;


% Colorize each axis for easy identification.
set(hax, 'XColor', [0.8 0 0.5], 'YColor', [0.1 0.7 0.3], 'ZColor', [0 0.4 0.7]);


% Set up the legend (unless there are too many clusters, in which case just show noise).
clustnames = cellstr(num2str(show));
if (exist('h_out', 'var')),  h_pts = [h_out, h_pts];   clustnames = cat(1, {'Outliers'}, clustnames);  end;

%if (numclusts < 12  &&  numclusts > 1)
if (length(show) < 18  &&  length(show) > 1)  
    hleg = legend(h_pts, clustnames, 0);    
elseif (exist('h_out', 'var')),
    hleg = legend(h_pts(1), clustnames{1}, 0);
end

%set(hfig, 'ButtonDownFcn', {@make_density, hax});
figure(hfig);  % bring figure to top



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

