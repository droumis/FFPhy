

%{

compatible with multitetrode anal
- multitetrode anal feeds all data in as varargins

TODO:
- add consensus ripple filtered LFP band
- add single unit spikes colored, organized, by subarea

LATER
- add multiunit spikes, color, organized by subarea
- add clusterless decode
- add clustered decode
%}

function out = dfa_plotDataChunks(idx, excludeIntervals, varargin)
appendindex = 1;
eventtype = 'rippleskons';
savefigs = 1;
pausefigs = 0;
pconf = paramconfig;
minstdthresh = 3;
maxvelocity = 4;
if ~isempty(varargin)
    assign(varargin{:});
end

day = idx(1,1);
epoch = idx(1,2);
% get the ripplekons events, which could be variably named
evid = find(contains(varargin(1:2:end), 'rippleskon'));
o = [1:2:length(varargin)]+1;
events = varargin{o(evid)};

evid = find(contains(varargin(1:2:end), 'noisekon'));
o = [1:2:length(varargin)]+1;
noise = varargin{o(evid)};

% out.excludeperiods = excludeperiods;
out.index = idx;
out.data = [];
out.dims = {};

% check for events
eventTime = [];
try
    eventTime = [events{day}{epoch}{1}.starttime events{day}{epoch}{1}.endtime];
catch
    fprintf('no events detected for day%d ep%d\n', day,epoch)
    return
end
% print proportion of included events
ecbefore = size(eventTime,1);
eventTime = eventTime(~isExcluded(eventTime(:,1),excludeIntervals),:);
ecafter = size(eventTime,1);
fprintf('%d / %d events included : d%d e%d\n', ecafter, ecbefore, day,epoch)
if isempty(eventTime)
    return
end
%% data intervals
validntrodes = find(cell2mat(cellfun(@(x) ~isempty(x), eeg{day}{epoch}, 'un', 0)));
% get start end of epoch from lfp starttime, endtime
% create an array of start, end intervals
startTime = eeg{day}{epoch}{validntrodes(1)}.starttime;
endTime = eeg{day}{epoch}{validntrodes(1)}.endtime;
timeSplits = [startTime:splitSize:endTime]';
timeSplits = [timeSplits(1:end-1) timeSplits(2:end)]; % start end
% lfp time
srate = eeg{day}{epoch}{validntrodes(1)}.samprate;
nsamps = length(eeg{day}{epoch}{validntrodes(1)}.data);
lfptime = [(0:nsamps-1)/srate+startTime]'; %repmat([,1,length(ntrodes));
% pos time
timecol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'time'));
postime = pos{day}{epoch}.data(:,timecol);
linpostime = linpos{day}{epoch}.statematrix.time;

%% rip zscpew
zrippwr = ((events{day}{epoch}{1}.powertrace-nanmean(events{day}{epoch}{1}.powertrace))...
    ./nanstd(events{day}{epoch}{1}.powertrace))';
zNoiseRipRatio = ((noise{day}{epoch}{1}.powertrace-nanmean(noise{day}{epoch}{1}.powertrace))...
    ./nanstd(noise{day}{epoch}{1}.powertrace))';

%% get nt's for each area
validntrodes = find(cell2mat(cellfun(@(x) ~isempty(x), eeg{day}{epoch}, 'un', 0)));
andef = animaldef(animal);
ntinfo = loaddatastruct(andef{2}, andef{3}, 'tetinfo');
%     reftet_keys = evaluatefilter(ntinfo, 'isequal($area, ''ref'')');
ca1Nts = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''ca1'')');
mecSpNts = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''mec'') && (isequal($subarea, ''supf'') || isequal($subarea, ''nsupf''))');
mecDpNts = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''mec'') && (isequal($subarea, ''deep'') || isequal($subarea, ''ndeep''))');

%% ---- reorder the nt's based on area,subarea ----
areaTags = cellfun(@(x) ntinfo{day}{epoch}{x}.area,num2cell(validntrodes),'un',0)';
subareaTags = cellfun(@(x) ntinfo{day}{epoch}{x}.subarea,num2cell(validntrodes),'un',0)';
ordernts = {'ca1', 'd'; ...
    'mec', 'deep'; 'mec', 'ndeep'; ...
    'mec', 'supf'; 'mec', 'nsupf'};
areaTagsNT = {};
subareaTagsNT = {};
ntsort = []; areabounds = []; areaid = {};
for i = 1:size(ordernts,1)
    ntinArea = find(all([ismember(areaTags, ordernts{i,1}) ismember(subareaTags, ...
        ordernts{i,2})],2));
    ntsort = [ntsort; ntinArea];
    areaid = [areaid; repmat(ordernts(i,:), length(ntinArea),1)];
end
ca1s = cellfun(@(x) any(strcmp(x,'d')), areaid(:,2),'un', 1);
dmecs = cellfun(@(x) any(strfind(x,'deep')), areaid(:,2),'un', 1);
smecs = cellfun(@(x) any(strfind(x,'supf')), areaid(:,2),'un', 1);
areaids = sum([1.*ca1s 2*dmecs 3*smecs],2);
ntsort = flipud(ntsort);
icolors = {};
icolors = colorPicker(areaTags, subareaTags);
icolors = icolors(ntsort);
%% GET SU
su_f = 'all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''),$tags,''un'',0))))';
su_keys = evaluatefilter(cellinfo{day}{epoch}, su_f);
su_keys = su_keys(ismember(su_keys(:,1), validntrodes),:);
if ~isempty(su_keys)
    [~,s] = ismember(su_keys(:,1), validntrodes(ntsort));
    [r, ri] = sort(s);
    sortsukeys = su_keys(ri,:);
    % create array (ntrode, cluster, spiketime, cellnum)
    supspikes = cellfun(@(x,y) [...
        repmat(x(1),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
        repmat(x(2),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
        spikes{day}{epoch}{x(1)}{x(2)}.data ...
        repmat(y,length(spikes{day}{epoch}{x(1)}{x(2)}.data),1)], ...
        num2cell(sortsukeys,2), ...
        num2cell(1:size(sortsukeys,1))', 'un', 0);
    sustack = cell2mat(supspikes(find(cell2mat(cellfun(@(x) ~isempty(x), supspikes, 'un', 0)))));
    sustack = sortrows(sustack,3);
else
    sustack = [];
end
%% MU
mu_keys = evaluatefilter(cellinfo{day}{epoch}, 'isequal($tags,{''mua''})');
% create array (ntrode, cluster, spiketime,ntsorted)
mupspikes = cellfun(@(x) [...
    repmat(x(1),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
    repmat(x(2),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
    spikes{day}{epoch}{x(1)}{x(2)}.data ...
    repmat(find(validntrodes(ntsort)==x(1)),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1)], ...
    num2cell(mu_keys,2), 'un', 0);
mustack = cell2mat(mupspikes(find(cell2mat(cellfun(@(x) ~isempty(x)&&size(x,2)==4, mupspikes, 'un', 0)))));
mustack = sortrows(mustack,3);
%% PLOT
Pp=load_plotting_params({'defaults','dataExplore'});
nrows = Pp.nrows;
for ts = 1:length(timeSplits(:,1))
    tstart = timeSplits(ts,1);
    tend = timeSplits(ts,2);
    % get ripple times in window
    swrInWinIdx = all([eventTime(:,1)>tstart eventTime(:,2)<tend],2);
    if ~any(swrInWinIdx)
        fprintf(sprintf('no rips in win %1.f : %1.f\n', tstart, tend)); continue; end
    if savefigs && ~pausefigs; close all;
        ifig = figure('Visible','off','units','normalized','position', Pp.position);
    else; ifig = figure('units','normalized','position',Pp.position); end
    set(gcf,'color','white')
    swrInWin = eventTime(swrInWinIdx,:);
    xs = swrInWin(:,1);
    xe = swrInWin(:,2);
    %% plot rip consensus power trace
    rippwrTime = events{day}{epoch}{1}.eegtimesvec_ref';
    pwrIdx = knnsearch(rippwrTime,[tstart; tend]);
    rippwrinWin = zrippwr(pwrIdx(1):pwrIdx(2));
    noiseinWin = zNoiseRipRatio(pwrIdx(1):pwrIdx(2));
    sfl = subaxis(nrows,1,1, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
        Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
    plot(rippwrTime(pwrIdx(1):pwrIdx(2))', noiseinWin,'LineWidth',2,'Color',[.8 .3 .3 .5]);
    hold on;
    plot(rippwrTime(pwrIdx(1):pwrIdx(2))', rippwrinWin,'LineWidth',2,'Color',[.3 .3 .3]);
    h=zoom; h.Motion = 'horizontal'; h.Enable = 'on';
    p = pan; p.Motion = 'horizontal'; set(gca, 'YGrid', 'off', 'XGrid', 'on');
    ylabel('rippwr'); set(gca, 'xticklabel', []); set(gca,'TickDir','out');
    axis tight; yl = ylim; xl = xlim;
    patch([xl(1) xl(2) xl(2) xl(1)]', [yl(1) yl(1)  minstdthresh  minstdthresh]','k', ...
        'linestyle', ':', 'edgecolor', 'none', 'facealpha', Pp.threshFAlpha);
    patch([xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,length(xe)),'y', ...
        'FaceAlpha',Pp.ripFAlpha, 'edgecolor','none'); hold off
    %% Plot LFP
    sfl = subaxis(nrows,1,2:4, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
        Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
    
    lfpIdx = knnsearch(lfptime,[tstart; tend]);
    eegstack = cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2))',...
        eeg{day}{epoch}(validntrodes),'un',0)');
    gridmat = Yoffset*repmat((1:length(eegstack(:,1)))', 1,length(eegstack(1,:)));
    eeggrid = (eegstack(ntsort,:) + gridmat)';
    ntrodeOffsets = gridmat(:,1);
    plfp = plot(lfptime(lfpIdx(1):lfpIdx(2))', eeggrid,'LineWidth',1);
    h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p = pan; p.Motion = 'horizontal';
    ylabel('LFP'); hold on; axis tight; yl = ylim; set(gca, 'xticklabel', []);
    set(gca,'TickDir','out'); set(plfp, {'color'}, icolors);
    set(gca, 'YGrid', 'off', 'XGrid', 'on');
    
    a_sA = cellfun(@(x,y,z) sprintf('%s %s %d',x,y,z),areaTags(ntsort), ...
        subareaTags(ntsort),num2cell(validntrodes(ntsort)'),'un', 0);
    yticks([ntrodeOffsets]); yticklabels(a_sA);
    patch([xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,length(xe)),'y', ...
        'FaceAlpha',.15, 'edgecolor','none'); hold off
    
    %% Plot ripple power lfp
    sfl = subaxis(nrows,1,5:6, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
        'MarginBottom', Pp.MgBm);
    % gath envelope magnitude of riplfp from this window and Vert Offset for plotting
    eegstack = double(cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2),3),...
        ripple{day}{epoch}(validntrodes),'un',0)))';
    im = imagesc(lfptime(lfpIdx(1):lfpIdx(2))', [1:length(eegstack(:,1))]-.5, ...
        flipud(eegstack(ntsort,:)));
    colormap(sfl,Pp.rippwrcmap);
    h=zoom; h.Motion ='horizontal'; h.Enable = 'on'; p = pan; p.Motion = 'horizontal';
    ylabel('ripLFP'); hold on; axis tight; yl = ylim; set(gca, 'xticklabel', []);
    set(gca,'TickDir','out'); yticks([ntrodeOffsets]); yticklabels(a_sA);
    patch([xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,length(xe)),'y', ...
        'FaceAlpha',0, 'edgecolor','w', 'edgealpha', .5);
    ab = find(diff(areaids));
    line(xl, [ab'; ab'], 'Color',[1 1 1 Pp.areaSepAlpha],'LineWidth', Pp.areaSepWidth);
    line([xl(1)-.01 xl(1)-.01], [yl(2) ab(2)], 'color',colorPicker({'mec'},{'supf'}), ...
        'linewidth', 4);
    line([xl(1)-.01 xl(1)-.01], [ab(2) ab(1)], 'color',colorPicker({'mec'}, {'deep'}),...
        'linewidth', 4);
    hold off;
    
    %% Plot theta power lfp
    sfl = subaxis(nrows,1,7:8, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
        'MarginBottom', Pp.MgBm);
    % envelope magnitude of thetalfp in window with vertical offset
    try
    eegstack = double(cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2),3),...
        theta{day}{epoch}(validntrodes),'un',0)))';
    catch
        warning('no theta eeg for %d %d', day, epoch);
        return
    end
    im = imagesc(lfptime(lfpIdx(1):lfpIdx(2))', [1:length(eegstack(:,1))]-.5,...
        flipud(eegstack(ntsort,:)));
    %         c = colormap(flipud(hsv));
    colormap(sfl,Pp.thetapwrcmap);
    h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p = pan; p.Motion = 'horizontal';
    ylabel('thetaLFP'); hold on; axis tight; yl = ylim;
    set(gca, 'xticklabel', []);
    set(gca,'TickDir','out');
    yticks([ntrodeOffsets]); yticklabels(a_sA);
    patch([xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,length(xe)),'y', ...
        'FaceAlpha',0, 'edgecolor','k', 'edgealpha', .5);
    line(xl, [ab'; ab'], 'Color',[1 1 1 Pp.areaSepAlpha],'LineWidth', Pp.areaSepWidth);
    line([xl(1)-.01 xl(1)-.01], [yl(2) ab(2)], 'color',colorPicker({'mec'},{'supf'}), ...
        'linewidth', 4);
    line([xl(1)-.01 xl(1)-.01], [ab(2) ab(1)], 'color',colorPicker({'mec'}, {'deep'}),...
        'linewidth', 4);
    hold off;
    
    %% Plot theta phase 
    sfl = subaxis(nrows,1,9:10, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
        'MarginBottom', Pp.MgBm);
    % envelope magnitude of thetalfp in window with vertical offset
    eegstack = double(cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2),2),...
        theta{day}{epoch}(validntrodes),'un',0)))';
    %         eegstack = diff([eegstack(:,1) eegstack],[],2);
    im = imagesc(lfptime(lfpIdx(1):lfpIdx(2))', [1:length(eegstack(:,1))]-.5,...
        flipud(eegstack(ntsort,:)));
    %         c = colormap(flipud(hsv));
    colormap(sfl,pink);
    h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p = pan; p.Motion = 'horizontal';
    ylabel('thetaLFP'); hold on; axis tight; yl = ylim;
    set(gca, 'xticklabel', []);
    set(gca,'TickDir','out');
    yticks([ntrodeOffsets]); yticklabels(a_sA);
    patch([xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,length(xe)),'y', ...
        'FaceAlpha',0, 'edgecolor','k', 'edgealpha', .5);
    line(xl, [ab'; ab'], 'Color',[1 1 1 Pp.areaSepAlpha],'LineWidth', Pp.areaSepWidth);
    line([xl(1)-.01 xl(1)-.01], [yl(2) ab(2)], 'color',colorPicker({'mec'},{'supf'}), ...
        'linewidth', 4);
    line([xl(1)-.01 xl(1)-.01], [ab(2) ab(1)], 'color',colorPicker({'mec'}, {'deep'}),...
        'linewidth', 4);
    hold off;
    
    %% Plot SU spikes
    sfp = subaxis(nrows,1,11:12,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
        Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
    if ~isempty(sustack)
        suinchunk = sustack(all([sustack(:,3)>tstart sustack(:,3)<tend],2),:);
        suclr = zeros(size(suinchunk,1),3);
        ca1spikes = ismember(suinchunk(:,1), ca1Nts);
        suclr(ca1spikes, :) = repmat(colorPicker({'ca1'}, {'d'}),sum(ca1spikes),1);
        deepmecspikes = ismember(suinchunk(:,1),mecDpNts);
        suclr(deepmecspikes, :) = repmat(colorPicker({'mec'}, {'deep'}),sum(deepmecspikes),1);
        supfmecspikes = ismember(suinchunk(:,1),mecSpNts);
        suclr(supfmecspikes, :) = repmat(colorPicker({'mec'}, {'supf'}),sum(supfmecspikes),1);
        a = scatter(suinchunk(:,3), suinchunk(:,4),Pp.SUsize,suclr,'filled');
        a.Marker = 'd';
         hold on; axis tight; yl = ylim; ylim([yl(1)-2 yl(2)+2]);
        yl = ylim; set(gca, 'xticklabel', []);
        set(gca,'TickDir','out'); set(gca, 'YGrid', 'off', 'XGrid', 'on');
        patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
            'FaceAlpha',.15, 'edgecolor','none'); hold off
    end
    ylabel('SU spikes');
    %% Plot MU spikes
    sfp = subaxis(nrows,1,13,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
        Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
    muinchunk = mustack(all([mustack(:,3)>tstart mustack(:,3)<tend],2),:);
    muclr = zeros(size(muinchunk,1),3);
    ca1spikes = ismember(muinchunk(:,1), ca1Nts);
    muclr(ca1spikes, :) = repmat(colorPicker({'ca1'}, {'d'}),sum(ca1spikes),1);
    deepmecspikes = ismember(muinchunk(:,1),mecDpNts);
    muclr(deepmecspikes, :) = repmat(colorPicker({'mec'}, {'deep'}),sum(deepmecspikes),1);
    supfmecspikes = ismember(muinchunk(:,1),mecSpNts);
    muclr(supfmecspikes, :) = repmat(colorPicker({'mec'}, {'supf'}),sum(supfmecspikes),1);
    a = scatter(muinchunk(:,3), muinchunk(:,4),1,muclr,'filled');
    a.Marker = 'd';
    ylabel('MU spikes'); hold on; axis tight; yl = ylim; ylim([yl(1)-1 yl(2)+1]);
    yl = ylim; set(gca, 'xticklabel', []);
    set(gca,'TickDir','out'); set(gca, 'YGrid', 'off', 'XGrid', 'on');
    patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
        'FaceAlpha',.15, 'edgecolor','none'); hold off
    %% plot velocity, pos, hd
    sfp = subaxis(nrows,1,14,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
        Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
    velcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'vel-loess'));
    xcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'x-loess'));
    ycol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'y-loess'));
    dircol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'dir-loess'));
    posIdx = knnsearch(postime(:,1),[tstart; tend]);
    vel = pos{day}{epoch}.data(posIdx(1):posIdx(2),velcol);
    veltime = postime(posIdx(1):posIdx(2));
    
    plot(veltime, vel, 'LineWidth', 2, 'Color', 'k', 'linestyle', '-');
    ylabel('2D Vel'); hold on; axis tight; set(gca, 'xticklabel', []);
    set(gca, 'YGrid', 'off', 'XGrid', 'on');
    set(gca,'TickDir','out'); yl = ylim; xl = xlim;
    velinds = vel <  maxvelocity;
    v = vel(velinds);
    v(v<4) = 1;
    v(v>=4) = 0;
    stem(veltime(velinds),v*yl(2),'.k','filled', 'linewidth', 5, 'color', [.8 .8 .8 .2]); hold on;
    plot(veltime, vel, 'LineWidth', 2, 'Color', 'k', 'linestyle', '-');
    %         line([xl(1) xl(2)],[4 4], 'linestyle', '--', 'color', 'b')
    patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
        'FaceAlpha',.15, 'edgecolor','none');
    hold off
    
    %% lin pos
    sfp = subaxis(nrows,1,15, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
        Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
    linposIdx = knnsearch(linpostime(:,1),[tstart; tend]);
    segmentnum = linpos{day}{epoch}.statematrix.segmentIndex(linposIdx(1):linposIdx(2));
    lindist = linpos{day}{epoch}.statematrix.linearDistanceToWells(...
        (linposIdx(1):linposIdx(2)),1);
    plot(linpostime(linposIdx(1):linposIdx(2)),lindist,'linewidth',3,'color',[.8 .8 .8]); %background pos
    xlabel('Time s'); ylabel('Lindist'); hold on; axis tight; ylim([0 200]); yl = ylim;
    set(gca,'TickDir','out'); set(gca, 'YGrid', 'on', 'XGrid', 'on');
    clr1 = [.6 .5 .2; .8 .5 .8; .6 0 .6; 0 0.4980 0.3255; 0 0.8 0.5216];
    for iseg = 1:5
        inds = []; inds = double(segmentnum == iseg); inds(find(inds==0))=nan;
        plot(linpostime(linposIdx(1):linposIdx(2)).*inds,lindist.*inds,'-', ...
            'linewidth',2,'color',clr1(iseg,:));
    end
    patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
        'FaceAlpha',.15, 'edgecolor','none'); hold off
    
    %% link x axis
    allAxesInFigure = findall(ifig,'type','axes'); linkaxes(allAxesInFigure, 'x');
    % ---- super axis -----
    sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
    sprtit = sprintf('DataChunk %s %d %d %1.f-%1.f', animal, day, epoch, tstart, tend);
    iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
    set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment','center','FontSize',12);
    % ---- pause, save figs ----
    if pausefigs; pause; end
    if savefigs; save_figure(sprintf('%s/dfa_plotDataChunk/%s/%s',pconf.andef{4},animal), sprtit); end
end