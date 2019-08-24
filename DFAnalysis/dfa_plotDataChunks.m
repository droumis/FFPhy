

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
savefigs = 0;
pausefigs = 1;
pconf = paramconfig;
if ~isempty(varargin)
    assign(varargin{:});
end

day = idx(1,1);
epoch = idx(1,2);
% get the ripplekons events, which could be variably named
evid = find(contains(varargin(1:2:end), 'rippleskon'));
o = [1:2:length(varargin)]+1;
events = varargin{o(evid)};

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

% get start end of epoch from lfp starttime, endtime
% create an array of start, end intervals
startTime = eeg{day}{epoch}{1}.starttime;
endTime = eeg{day}{epoch}{1}.endtime;
timeSplits = [startTime:splitSize:endTime]';
timeSplits = [timeSplits(1:end-1) timeSplits(2:end)]; % start end
% recreate lfpepoch time
srate = eeg{day}{epoch}{1}.samprate;
nsamps = length(eeg{day}{epoch}{1}.data);
ntrodes = find(cell2mat(cellfun(@(x) ~isempty(x), eeg{day}{epoch}, 'un', 0)));
lfptime = [(0:nsamps-1)/srate+startTime]'; %repmat([,1,length(ntrodes));
timecol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'time'));
%postime
postime = pos{day}{epoch}.data(:,timecol);
linpostime = linpos{day}{epoch}.statematrix.time;

%get su idx

data_keys = cellfetch(cellinfo{day}{epoch}, '');
data_keys = data_keys.index;
mu_keys = evaluatefilter(cellinfo{day}{epoch}, 'isequal($tags, {''mua''})');
su_keys = data_keys(~ismember(data_keys, mu_keys, 'rows'),:);
supspikes = cellfun(@(x,y) [repmat(x(1),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
    repmat(x(2),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
    spikes{day}{epoch}{x(1)}{x(2)}.data ...
    repmat(y,length(spikes{day}{epoch}{x(1)}{x(2)}.data),1)], num2cell(su_keys,2),num2cell(1:size(su_keys,1))', 'un', 0);
sustack = cell2mat(supspikes(find(cell2mat(cellfun(@(x) ~isempty(x), supspikes, 'un', 0)))));
sustack = sortrows(sustack,3);
% get nt's for each area
andef = animaldef(animal);
ntinfo = loaddatastruct(andef{2}, andef{3}, 'tetinfo');
reftet_keys = evaluatefilter(ntinfo, 'isequal($area, ''ref'')');
ca1tet_keys = evaluatefilter(ntinfo, 'isequal($area, ''ca1'')');
adet_mecSupfNts = evaluatefilter(ntinfo, 'isequal($area, ''mec'') && isequal($subarea, ''supf'')');
adet_mecDeepNts = evaluatefilter(ntinfo, 'isequal($area, ''mec'') && isequal($subarea, ''deep'')');
refNts = reftet_keys(ismember(reftet_keys(:,[1 2]),[day epoch], 'rows'),3);
ca1Nts = ca1tet_keys(ismember(ca1tet_keys(:,[1 2]), [day epoch], 'rows'),3);
mecSpNts = adet_mecSupfNts(ismember(adet_mecSupfNts(:,[1 2]), [day epoch], 'rows'),3);
mecDpNts = adet_mecDeepNts(ismember(adet_mecDeepNts(:,[1 2]), [day epoch], 'rows'),3);

%% ---- reorder the LFP traces by suparea|area|subarea tags (in that priority) ----

ntinfoAllTags = cellfetch(ntinfo, '', 'alltags', 1);
ntrodesIdx = ntinfoAllTags.index(ismember(ntinfoAllTags.index(:,[1 2]), [day epoch], 'rows'),:);
[~, tagIndMap] = ismember(ntrodesIdx,ntinfoAllTags.index, 'rows');
ntrodeTags = ntinfoAllTags.values(tagIndMap);
try
    numsumSupAreas = cellfun(@(x) sum(uint16(x.suparea)), ntrodeTags, 'un', 0);
    numsumAreas = cellfun(@(x) sum(uint16(x.area)), ntrodeTags, 'un', 0);
    numsumSubAreas = cellfun(@(x) sum(uint16(x.subarea)), ntrodeTags, 'un', 0);
    strSupAreas = cellfun(@(x) x.suparea, ntrodeTags, 'un', 0);
    strAreas = cellfun(@(x) x.area, ntrodeTags, 'un', 0);
    strSubAreas = cellfun(@(x) x.subarea, ntrodeTags, 'un', 0);
catch
    error('all ntrodes need to have a suparea, subarea, and area tag, even if blank')
end
icolors = colorPicker(strAreas, strSubAreas);
numsumallareatags = cell2mat([numsumSupAreas numsumAreas numsumSubAreas]);
[numsumallSort, numsumSortInds] = sortrows(numsumallareatags);%,[-1 -2 -3]); % -Col to 'descend'
icolors = icolors(numsumSortInds,:);

for ts = 1:length(timeSplits(:,1))
    tstart = timeSplits(ts,1);
    tend = timeSplits(ts,2);
    % get ripple times in window
    swrInWinIdx = all([eventTime(:,1)>tstart eventTime(:,2)<tend],2);
    if ~any(swrInWinIdx)
        fprintf(sprintf('no rips in win %1.f : %1.f\n', tstart, tend));
        continue
    end
    swrInWin = eventTime(swrInWinIdx,:);
    xs = swrInWin(:,1);
    xe = swrInWin(:,2);
    
    %% Plot LFP
    Pp=load_plotting_params({'defaults','dataExplore'});
    % ---- init fig----
    if savefigs && ~pausefigs; close all; 
        ifig = figure('Visible','off','units','normalized','position', Pp.position);
    else; ifig = figure('units','normalized','position',Pp.position); end 
    set(gcf,'color','white')
    % plot LFP
    % get time split indices into lfptime
    lfpIdx = knnsearch(lfptime,[tstart; tend]);
    % gath lfp from this window and Vert Offset for plotting
    eegstack=cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2))',eeg{day}{epoch},'un',0)');
    a = eegstack(numsumSortInds,:);
    gridmat = Yoffset*repmat((1:length(eegstack(:,1)))', 1,length(eegstack(1,:)));
    eeggrid = (a + gridmat)';
    ntrodeOffsets = gridmat(:,1);
    sfl = subaxis(6,1,[1:3], 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
        Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
    hold on;
    plfp = plot(lfptime(lfpIdx(1):lfpIdx(2))', eeggrid,'LineWidth',1,'Color',[.2 .2 .2]);
    h = zoom; h.Motion = 'horizontal'; h.Enable = 'on'; axis tight; 
    set(gca, 'xtick', []); ylabel('lfp'); p = pan; p.Motion = 'horizontal'; yl = ylim;
    % set the colors of the lfp traces by the area, subarea
    set(plfp, {'color'}, icolors); %num2cell(lfpC, 2));
    a_sA = cellfun(@(x,y,z) sprintf('%s %s %d',x,y,z),strAreas(numsumSortInds), strSubAreas(numsumSortInds),num2cell(numsumSortInds),'un', 0);
    yticks([ntrodeOffsets]); yticklabels(a_sA);
    patch([xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,length(xe)),'y', ...
        'FaceAlpha',.15, 'edgecolor','none');
    hold off
    
    %% Plot SU spikes
    sfp = subaxis(6,1,4,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
        Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
    suinchunk = sustack(all([sustack(:,3)>tstart sustack(:,3)<tend],2),:);
    scatter(suinchunk(:,3), suinchunk(:,4),10,sz,'filled');
    axis tight; set(gca, 'xtick', [])
    
    %% plot velocity, pos, hd
    velcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'vel-loess'));
    xcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'x-loess'));
    ycol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'y-loess'));
    dircol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'dir-loess'));
    % get time split indices into pos
    posIdx = knnsearch(postime(:,1),[tstart; tend]);
    % velocity
    sfp = subaxis(6,1,5,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
        Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
    plot(postime(posIdx(1):posIdx(2))', pos{day}{epoch}.data(posIdx(1):posIdx(2),velcol), ...
        'LineWidth', 2, 'Color', [.1 .1 .1])
    axis tight; set(gca, 'xtick', []); ylabel('2D Velocity'); yl = ylim; hold off
    patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
        'FaceAlpha',.15, 'edgecolor','none'); hold on;
%     
%     % angular velocity
%     sfp = subaxis(6,1,5,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
%         Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
%     plot(postime(posIdx(1)+1:posIdx(2))', ...
%         abs(diff(pos{day}{epoch}.data(posIdx(1):posIdx(2),dircol))), ...
%         'LineWidth', 2, 'Color', [.5 .5 .5], 'LineStyle', ':', 'MarkerSize', .5)
%     axis tight; set(gca, 'xtick', []); xlabel('time s'); yl = ylim; hold off
%     patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
%         'FaceAlpha',.15, 'edgecolor','none'); hold on;
%     
%     % lin pos
    sfp = subaxis(6,1,6, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
        Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
%     plot(postime(posIdx(1):posIdx(2))', pos{day}{epoch}.data(posIdx(1):posIdx(2),[xcol ycol]), ...
%         'LineWidth', 2)
    linposIdx = knnsearch(linpostime(:,1),[tstart; tend]);
    segmentnum = linpos{day}{epoch}.statematrix.segmentIndex(linposIdx(1):linposIdx(2));
    lindist = linpos{day}{epoch}.statematrix.linearDistanceToWells((linposIdx(1):linposIdx(2)),1);
    plot(linpostime(linposIdx(1):linposIdx(2)),lindist,'linewidth',2,'color',[.8 .8 .8]); %background pos
    hold on; axis tight; xlabel('time s'); ylabel('lindist'); ylim([0 200]);yl = ylim;
    clr1 = [.6 .5 .2; .8 .5 .8; .6 0 .6; 0 0.4980 0.3255; 0 0.8 0.5216];
    %             line(postime, lindist, segmentnum);
    for iseg = 1:5
        inds = [];
        inds = double(segmentnum == iseg);
        inds(find(inds==0))=nan;
        plot(linpostime(linposIdx(1):linposIdx(2)).*inds,lindist.*inds,'-','linewidth',2,'color',clr1(iseg,:));
        hold on;
    end 
    patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
        'FaceAlpha',.15, 'edgecolor','none'); hold on;
    hold off
    % link x axis
    allAxesInFigure = findall(ifig,'type','axes');
    linkaxes(allAxesInFigure, 'x');
    %     xlim([startTime, startTime+20])
    
    % ---- super axis -----
    sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
    sprtit = sprintf('DataChunk %s %d %d %1.f-%1.f', animal, day, epoch, tstart, tend);
    iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
    set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment','center','FontSize',12);
    % ---- pause, save figs ----
    if pausefigs; pause; end
    if savefigs; save_figure([pconf.andef{4},'/dfa_plotDataChunk/'],sprtit); end
end
% end

