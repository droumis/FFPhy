

function out = dfa_lickXCorrSpikes(index, excludeperiods, varargin)
fprintf('%d %d %d %d\n',index)
reqData = {'lick', 'spikes', 'cellinfo'};
for s = 1:length(reqData)
    if ~any(cellfun(@(x) strcmp(x,reqData{s}), varargin(1:2:end), 'un', 1))
        error('missing required data')
    end
end

pconf = paramconfig;
bin = .02;
tmax = 1;
eventName = 'lick';
psthbin = .001;
shuff = 100; % ms
TF = 1; % legacy tetfilter
plotfigs = 1;
displayplots = 0;
saveplots = 1;
savefigas = {'png', 'mfig'};
if ~isempty(varargin)
    assign(varargin{:});
end

% init out
out.index = index;
out.dims = {};
out.data = [];

day = index(1);
epoch = index(2);
ntrode = index(3);
cluster = index(4);
try
    spiketimes = spikes{day}{epoch}{ntrode}{cluster}.data(:,1);
catch
    return
end

% get events
evvar = varargin{find(cellfun(@(x) ~isempty(x), strfind(varargin(1:2:end), ...
    eventName), 'un', 1))*2-1};
events = eval(evvar);
try
    eventTimes = events{day}{epoch}.starttime;
    id = events{day}{epoch}.id;
catch
    eventTimes = events{day}{epoch}{TF}.starttime;
    id = events{day}{epoch}{TF}.id;
end
% timefilter events
includetimes = ~isExcluded(eventTimes, excludeperiods); % result: 1 include, 0 exclude
eventTimes = eventTimes(includetimes);
% if isempty(eventtimes)
%     return
% end
realxcs = spikexcorr(spiketimes,eventTimes, bin, tmax);
% time shift shuffle
x = [];
shiftby = randi([-shuff shuff],1,1000)/1e3; %ms
for i = 1:1000
    %         Y = circshift(licks(:,1),shiftby);
    lickShift = sort(eventTimes+shiftby(i));
    xcs = spikexcorr(spiketimes, lickShift, bin, tmax);
    x(i,:) = xcs.c1vsc2;
end
xstd = std(x,[],1);
xmean = mean(x,1);
if plotfigs
    %% plot
    Pp=load_plotting_params({'defaults','dfa_lickXCorrSpikes'});
    if saveplots && ~displayplots; close all;
        ifig = figure('Visible','off','units','normalized','position', Pp.position, ...
            'color','white', 'InvertHardcopy', 'off');
    else
        ifig = figure('units','normalized','position',Pp.position, 'color','white', ...
            'InvertHardcopy', 'off');
    end
    
    nrow = 3;
    ncol = 1;
    %% lick psth
    psthtime = (-tmax-0.5*psthbin):psthbin:(tmax+0.5*psthbin);
    if ~isempty(eventTimes)
        psth = nan(length(eventTimes),length(psthtime));
        for r=1:length(eventTimes)
            shist = histc(spiketimes, eventTimes(r) + psthtime);
            psth(r,:) = logical(shist);
        end
    end
    sf1 = subaxis(nrow,ncol,1,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
        'MarginBottom', Pp.MgBm); set(gca, 'Tag', 'ripkons');
    h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p=pan; p.Motion='horizontal';
    [xx, yy] = find(psth');
    uclr = cell2mat(cellfun(@(x) repmat(x(1),x(2),1),num2cell([id sum(psth,2)],2), ...
        'un',0));
    f = scatter(xx/1000-1.001,yy,Pp.psthSize, uclr, 'filled');
    colormap(cbrewer('qual', 'Accent', 3))
    set(gca, 'color', 'k')
    sf1.GridColor = 'w';
    f.Marker = 'd';
    f.MarkerEdgeAlpha = 0;
    f.MarkerFaceAlpha = 0.4;
    set(gca, 'TickDir','out', 'YGrid', 'off', 'XGrid', 'on', 'TickLength', [0.001 0], ...
        'Tag', 'ndata', 'Xticklabel',[]);
    axis tight;
    ylabel('licknum','FontSize',12,'FontWeight','bold', 'FontName','Arial')
    axis tight
    line([0 0],ylim, 'color','red', 'linewidth',2, 'Color', [1 0 0 .5]);
    %% psth pdf
    sf2 = subaxis(nrow,ncol,2);
    h = histogram(xx/1000-1.001, psthtime(1:20:end), 'facecolor', 'k', 'Normalization', ...
        'pdf'); 
    h.EdgeColor = 'none';
    axis tight
    line([0 0],ylim, 'color','red', 'linewidth',2, 'Color', [1 0 0 .5]);
    ylabel('pdf','FontSize',12,'FontWeight','bold', 'FontName','Arial')
    set(gca, 'YGrid', 'off', 'XGrid', 'on','TickDir','out', 'TickLength', [0.001 0]);
    xticklabels([])
    hold off;
    %% lick XCORR spikes
    sf3 = subaxis(nrow,ncol,3);
    realxcs = spikexcorr(spiketimes,eventTimes, bin, tmax);
    plot(realxcs.time, realxcs.c1vsc2', 'k', 'linewidth', 2)%, 'facecolor', 'k')%, '-k', 'filled', 'markersize', 5); 
    hold on;
    axis tight;
    plot(xcs.time, xmean, 'color', 'b', 'linewidth', 2); hold on;
    fill([xcs.time'; flipud(xcs.time')],[xmean'-xstd';flipud(xmean'+xstd')],'b',...
        'linestyle','none', 'facealpha', .1);
    line([0 0],ylim, 'color','red', 'linewidth',2, 'Color', [1 0 0 .5]);
    hold off
    ylabel('xcorr')
    xlabel('lag time s')
    set(gca, 'YGrid', 'off', 'XGrid', 'on','TickDir','out', 'TickLength', [0.001 0]);

    %% ----- link x axis -----
    allAxesInFigure = findall(ifig,'type','axes'); linkaxes(allAxesInFigure, 'x');
    % ---- super axis -----
    sprtit = sprintf('%s %d %d %d %d lickXCorrSpikes bin%.0f shuff%.0f',animal,day,epoch,...
        ntrode, cluster, bin*1e3, shuff);
    sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
    iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
    set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment','center', ...
        'FontSize',12);
    h = get(gcf,'Children');
    set(gcf,'Children',flip(h)); % super axis to bottom. allows for zoom/pan
    % ---- pause, save figs ----
    if displayplots; pause; end
    if saveplots
        [~, fname,~] = fileparts(mfilename('fullpath'));
        outdir = sprintf('%s/%s/%s/', pconf.andef{4},fname,animal);
        save_figure(outdir, sprtit, 'savefigas', savefigas); 
    end
end
end
