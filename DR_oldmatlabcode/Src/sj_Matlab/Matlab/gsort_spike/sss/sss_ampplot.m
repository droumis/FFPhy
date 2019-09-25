function sss_ampplot(spikes, assignments, show, nch, dopca)
% Amp plots AND PC PLOTS for diff channels of tetrode

% Shantanu
%   Last Modified By: sbm on Wed Aug 17 18:29:18 2005

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
%cmap = hot(2*length(show));
%end

if (nargin < 4), nch=4; end

%%numclusts = length(clusters);
%numclusts = length(show);   %%%% MY CHANGE: SHANTANU FEB 2006

show(show == 0) = [];

%%%% PRE-CALCULATE PCA FOR BIG DATASETS

% if dopca==1
%     % SVD the data because we want to use the PCs as axes.
%     [pca.scores,pca.u,pca.s,pca.v] = pcasvd(spikes.waveforms);
%     spikes.pca = pca;
%     data = pca.scores;
% end

if dopca==1
    data=spikes.pca.scores; %%% TAKE 1ST PC from sj_pcasvd
%     pca=spikes.pca;
%     data = pca.scores;
end

%%% TO NOT SHOW OUTLIERS
clusters(clusters == 0) = []; 
numclusts = length(clusters);


if nch==2

    for nch=1:2     %2 channels for stereotrode
        cmd=sprintf('ampl(:,nch)=abs(max(spikes.waveforms_ch%d,[],2)) + abs(min(spikes.waveforms_ch%d,[],2));',nch,nch); eval(cmd);
    end

    afig=figure; hold on;
    %    redimscreen; orient(gcf,'landscape');
    set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');
    ylabel('Electrode 2' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Electrode 1','FontSize', 14, 'FontWeight', 'bold');

    for clu = 1:numclusts
        sel = find(assignments == clusters(clu));
        inds{clu} = sel;

        hndl = plot(ampl(sel,1), ampl(sel,2), '.','MarkerSize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.1);
        elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
        else                                     set(hndl, 'Color', Clgy,'MarkerSize',0.1);
        end;
        title('AMP PLOTS')

        if (clusters(clu) == 0),                 h_out = hndl;
        elseif (ismember(clusters(clu), show)),  h_pts(clusters(clu)) = hndl;
        end
    end

    h_pts = h_pts(show);
    if (length(find(h_pts)) == 1), set(h_pts, 'MarkerSize', 2, 'Marker', '.');  end;
    uistack(h_pts, 'top');

    hold off;

    % Set up the legend (unless there are too many clusters, in which case just show noise).
    clustnames = cellstr(num2str(show));
    if (exist('h_out', 'var')),  h_pts = [h_out, h_pts];   clustnames = cat(1, {'Outliers'}, clustnames);  end;

    %if (numclusts < 12  &&  numclusts > 1)
    if (length(show) < 12  &&  length(show) > 1)    %% SHANTANU: FEB 2006: CHANGE
        hleg = legend(h_pts, clustnames, 0);
    elseif (exist('h_out', 'var')),
        hleg = legend(h_pts(1), clustnames{1}, 0);
    end

    set(afig, 'ButtonDownFcn', {@make_density, ampl});

    return
end

% Make Ampl_plot

for nch=1:4     %4 channels for tetrode
    cmd=sprintf('ampl(:,nch)=abs(max(spikes.waveforms_ch%d,[],2)) + abs(min(spikes.waveforms_ch%d,[],2));',nch,nch); eval(cmd);
end

% Make the plot.
figure;  hfig = gcf;  hax = gca;    hold on; redimscreen100(100);
%    redimscreen; orient(gcf,'landscape');
set(gcf, 'PaperPositionMode', 'auto');
set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');

%whitebg(hfig);

for clu = 1:numclusts
    sel = find(assignments == clusters(clu));
    inds{clu} = sel;

    subplot(3,2,1); hold on;
    hndl = plot(ampl(sel,1), ampl(sel,2), '.','MarkerSize',2);
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.1);
    elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
    else                                     set(hndl, 'Color', Clgy,'MarkerSize',0.01);
    end;
    title('AMP PLOTS','FontSize', 14, 'FontWeight', 'bold');
    ylabel('Electrode 2' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Electrode 1','FontSize', 14, 'FontWeight', 'bold');

    subplot(3,2,2); hold on;
    ylabel('Electrode 3' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Electrode 1','FontSize', 14, 'FontWeight', 'bold');
    hndl = plot(ampl(sel,1), ampl(sel,3), '.','MarkerSize',2);
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.1);
    elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
    else                                     set(hndl, 'Color', Clgy,'MarkerSize',0.01);
    end;

    subplot(3,2,3); hold on;
    ylabel('Electrode 4' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Electrode 1','FontSize', 14, 'FontWeight', 'bold');
    hndl = plot(ampl(sel,1), ampl(sel,4), '.','MarkerSize',2);
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.1);
    elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
    else                                     set(hndl, 'Color', Clgy,'MarkerSize',0.01);
    end;

    subplot(3,2,4); hold on;
    ylabel('Electrode 3' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Electrode 2','FontSize', 14, 'FontWeight', 'bold');
    hndl = plot(ampl(sel,2), ampl(sel,3), '.','MarkerSize',2);
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.1);
    elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
    else                                     set(hndl, 'Color', Clgy,'MarkerSize',0.01);
    end;

    subplot(3,2,5); hold on;
    ylabel('Electrode 4' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Electrode 2','FontSize', 14, 'FontWeight', 'bold');
    hndl = plot(ampl(sel,2), ampl(sel,4), '.','MarkerSize',2);
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.1);
    elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
    else                                     set(hndl, 'Color', Clgy,'MarkerSize',0.01);
    end;

    subplot(3,2,6); hold on;
    ylabel('Electrode 4' ,'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Electrode 3','FontSize', 14, 'FontWeight', 'bold');
    hndl = plot(ampl(sel,3), ampl(sel,4), '.','MarkerSize',2);
    %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
    if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.1);
    elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
    else                                     set(hndl, 'Color', Clgy,'MarkerSize',0.01);
    end;

    if (clusters(clu) == 0),                 h_out = hndl;
    elseif (ismember(clusters(clu), show)),  h_pts(clusters(clu)) = hndl;
    end

end

h_pts = h_pts(show);
if (length(find(h_pts)) == 1), set(h_pts, 'MarkerSize', 2, 'Marker', '.');  end;

hold off;

% Set up the legend (unless there are too many clusters, in which case just show noise).
clustnames = cellstr(num2str(show));
if (exist('h_out', 'var')),  h_pts = [h_out, h_pts];   clustnames = cat(1, {'Outliers'}, clustnames);  end;

%if (numclusts < 12  &&  numclusts > 1)
if (length(show) < 12  &&  length(show) > 1)    %% SHANTANU: FEB 2006: CHANGE
    hleg = legend(h_pts, clustnames, 0);
elseif (exist('h_out', 'var')),
    hleg = legend(h_pts(1), clustnames{1}, 0);
end


if dopca==1
    
    %%%%%% PRE-STORED PCAS %%%%%%%%%%%%%%% 
    %data=spikes.pca.scores;
    data=spikes.pca.scores; %%% TAKE 1ST PC from sj_pcasvd
    
    % Make the PCA plot.
    figure;  hfig = gcf;     hax = gca;    hold on; redimscreen100(100);
    %whitebg(hfig);
    %    redimscreen; orient(gcf,'landscape');
    set(gcf, 'PaperPositionMode', 'auto');
    set(0,'defaultaxesfontsize',14);set(0,'defaultaxesfontweight','bold');

    for clu = 1:numclusts
        sel = find(assignments == clusters(clu));
        inds{clu} = sel;

        subplot(3,2,1); hold on;
        ylabel('PCA: Ch2' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('PCA: Ch1','FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(data(sel,1), data(sel,2), '.','MarkerSize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.5);
        elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
        else                                     set(hndl, 'Color', Clgy,'MarkerSize',1);
        end;
        title('PCA PLOTS')

        subplot(3,2,2); hold on;
        ylabel('PCA: E3' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('PCA: E1','FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(data(sel,1), data(sel,3), '.','MarkerSize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.5);
        elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
        else                                     set(hndl, 'Color', Clgy,'MarkerSize',1);
        end;

        subplot(3,2,3); hold on;
        ylabel('PCA: E4' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('PCA: E1','FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(data(sel,1), data(sel,4), '.','MarkerSize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.5);
        elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
        else                                     set(hndl, 'Color', Clgy,'MarkerSize',1);
        end;

        subplot(3,2,4); hold on;
        ylabel('PCA: E3' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('PCA: E2','FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(data(sel,2), data(sel,3), '.','MarkerSize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.5);
        elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
        else                                     set(hndl, 'Color', Clgy,'MarkerSize',1);
        end;

        subplot(3,2,5); hold on;
        ylabel('PCA: E4' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('PCA: E2','FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(data(sel,2), data(sel,4), '.','MarkerSize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.5);
        elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
        else                                     set(hndl, 'Color', Clgy,'MarkerSize',1);
        end;

        subplot(3,2,6); hold on;
        ylabel('PCA: E4' ,'FontSize', 14, 'FontWeight', 'bold');
        xlabel('PCA: E3','FontSize', 14, 'FontWeight', 'bold');
        hndl = plot(data(sel,3), data(sel,4), '.','MarkerSize',2);
        %set(hndl, 'ButtonDownFcn', {@raise_me, hndl});
        if (clusters(clu) == 0),                 set(hndl, 'Color', [0 0 0], 'Marker', 'x','MarkerSize',0.5);
        elseif (ismember(clusters(clu), show)),  set(hndl, 'Color', cmap(clu,:));
        else                                     set(hndl, 'Color', Clgy,'MarkerSize',2);
        end;

        if (clusters(clu) == 0),                 h_out = hndl;
        elseif (ismember(clusters(clu), show)),  h_pts(clusters(clu)) = hndl;
        end

    end

    h_pts = h_pts(show);
    if (length(find(h_pts)) == 1), set(h_pts, 'MarkerSize', 2, 'Marker', '.');  end;

    hold off;


    % Set up the legend (unless there are too many clusters, in which case just show noise).
    clustnames = cellstr(num2str(show));
    if (exist('h_out', 'var')),  h_pts = [h_out, h_pts];   clustnames = cat(1, {'Outliers'}, clustnames);  end;

    %if (numclusts < 12  &&  numclusts > 1)
    if (length(show) < 12  &&  length(show) > 1)    %% SHANTANU: FEB 2006: CHANGE
        hleg = legend(h_pts, clustnames, 0);
    elseif (exist('h_out', 'var')),
        hleg = legend(h_pts(1), clustnames{1}, 0);
    end
end


%%%%%%%%%%% Density Plot for Stereotrode Amplitude  %%%%%%%%%%%%%%%%

function make_density(afig,event,ampl)
figure; hdens = gca;
histxy(ampl(:,1),ampl(:,2),250,1);
%set(gca,properties2d{:},'Color',[0 0 0.508]);
HGlogmenu(findobj(gca,'Type','Image'));