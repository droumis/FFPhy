

%{

compatible with multitetrode anal
- multitetrode anal feeds all data in as varargins

LATER
- add clusterless decode
- add clustered decode
%}

function out = dfa_plotDataChunks(idx, excludeIntervals, varargin)
pconf = paramconfig;
savefigs = 1;
pausefigs = 0;
minstdthresh = 3;
maxvelocity = 4;
savefigas = 'png';
centerEvents = 0;
splitSize = 30;
Yoffset = 600;
skipNoRipChunks = 0;
skipNoEventChunks = 0;
if ~isempty(varargin)
    assign(varargin{:});
end

day = idx(1,1);
epoch = idx(1,2);
% get the ripplekons events, which could be variably named
evid = find(contains(varargin(1:2:end), 'rippleskon'));
o = [1:2:length(varargin)]+1;
events = varargin{o(evid)};

% evid = find(contains(varargin(1:2:end), 'noisekon'));
% o = [1:2:length(varargin)]+1;
% noise = varargin{o(evid)};

% out.excludeperiods = excludeperiods;
out.index = idx;
out.data = [];
out.dims = {};
try
    eventTime = [events{day}{epoch}{1}.starttime events{day}{epoch}{1}.endtime];
catch
    fprintf('no events detected for day%d ep%d\n', day,epoch)
    return
end
eventTime = eventTime(~isExcluded(eventTime(:,1),excludeIntervals),:);
% data intervals for tetrode already filtered for inclusion
validnt = find(cell2mat(cellfun(@(x) ~isempty(x), eeg{day}{epoch}, 'un', 0)));
% get start end of epoch from lfp starttime, endtime
% create an array of start, end intervals
startTime = eeg{day}{epoch}{validnt(1)}.starttime;
endTime = eeg{day}{epoch}{validnt(1)}.endtime;
if centerEvents
    timeSplits = [eventTime(:,1)-(splitSize/2) eventTime(:,2)+(splitSize/2)];
    % exclude chunks out of range
    while timeSplits(1,1) < startTime
        timeSplits(1,:) = [];
    end
    while timeSplits(end,2) > endTime
        timeSplits(end,:) = [];
    end
else
    if splitSize == 0
        timeSplits = [startTime endTime];
    else
        timeSplits = [startTime:splitSize:endTime]';
        timeSplits = [timeSplits(1:end-1) timeSplits(2:end)]; % start end
    end
end
% lfp time
srate = eeg{day}{epoch}{validnt(1)}.samprate;
nsamps = length(eeg{day}{epoch}{validnt(1)}.data);
lfptime = [(0:nsamps-1)/srate+startTime]'; %repmat([,1,length(ntrodes));

% pos time
timecol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'time'));
postime = pos{day}{epoch}.data(:,timecol);
linpostime = linpos{day}{epoch}.statematrix.time;

% stack lindist outer arm
lindist1D = linpos{day}{epoch}.statematrix.linearDistanceToWells(:,1);
outerIdx = ismember(linpos{day}{epoch}.statematrix.segmentIndex,[4,5]);
innermax = max(lindist1D(ismember(linpos{day}{epoch}.statematrix.segmentIndex,[2,3])));
centermax = max(lindist1D(ismember(linpos{day}{epoch}.statematrix.segmentIndex,1)));
if isempty(centermax)
    centermax = 0;
end
lindist1D(outerIdx) = lindist1D(outerIdx)+(innermax-centermax);
% get lindist segments start, end
segSE = []; c = 0;
for iseg = 1:length(unique(linpos{day}{epoch}.statematrix.segmentIndex))
    segstart = min(lindist1D(linpos{day}{epoch}.statematrix.segmentIndex == iseg));
    segend = max(lindist1D(linpos{day}{epoch}.statematrix.segmentIndex == iseg));
    if ~isempty(segstart)
        c = c+1;
        segSE(c,1) = segstart;
        segSE(c,2) = segend;
    end
end

% lfp ripple zcored pwr
zrippwr = ((events{day}{epoch}{1}.powertrace-nanmean(events{day}{epoch}{1}.powertrace))...
    ./nanstd(events{day}{epoch}{1}.powertrace))';
% 4:12Hz filtered ripple zcored pwr
bpzpwr = zscore(bpfilt(events{day}{epoch}{1}.powertrace,4, 12, 1500, 0));
%     bpzpwr = bandpass(zrippwr,[4 12],1500);
%     zNoiseRipRatio = ((noise{day}{epoch}{1}.powertrace-nanmean(noise{day}{epoch}{1}.powertrace))...
%         ./nanstd(noise{day}{epoch}{1}.powertrace))';

% process ntinfo ino nt's for each area 
andef = animaldef(animal);
ntinfo = loaddatastruct(andef{2}, andef{3}, 'tetinfo');
%     reftet_keys = evaluatefilter(ntinfo, 'isequal($area, ''ref'')');
c1Nt = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''ca1'')');
smNt = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''mec'') && (isequal($subarea, ''supf'') || isequal($subarea, ''nsupf''))');
dmNt = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''mec'') && (isequal($subarea, ''deep'') || isequal($subarea, ''ndeep''))');
% ---- reorder the nt's based on area,subarea ----
arTag = cellfun(@(x) ntinfo{day}{epoch}{x}.area,num2cell(validnt),'un',0)';
sbTag = cellfun(@(x) ntinfo{day}{epoch}{x}.subarea,num2cell(validnt),'un',0)';
ordernts = {'ca1', 'd'; ...
    'mec', 'deep'; 'mec', 'ndeep'; ...
    'mec', 'supf'; 'mec', 'nsupf'; ...
    'ref', 'ca1'; 'ref', 'mec'};
% determine the order of areas
ntsort = []; areaid = {};
for i = 1:size(ordernts,1)
    ntinArea = find(all([ismember(arTag, ordernts{i,1}) ismember(sbTag, ...
        ordernts{i,2})],2));
    ntsort = [ntsort; ntinArea];
    areaid = [areaid; repmat(ordernts(i,:), length(ntinArea),1)];
end

ca1s = cellfun(@(x) any(strcmp(x,'d')), areaid(:,2),'un', 1);
dmecs = cellfun(@(x) any(strfind(x,'deep')), areaid(:,2),'un', 1);
smecs = cellfun(@(x) any(strfind(x,'supf')), areaid(:,2),'un', 1);
areaids = sum([1.*ca1s 2*dmecs 3*smecs],2);
ntsort = flipud(ntsort);

% sorted colors
icolors = colorPicker(arTag, 'subtags', sbTag); icolors = icolors(ntsort);

% Process SU into sustack (ntrode, cluster, spiketime, cellnum)
su_f = 'all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''),$tags,''un'',0))))';
su_keys = evaluatefilter(cellinfo{day}{epoch}, su_f);
su_keys = su_keys(ismember(su_keys(:,1), validnt),:);
[~,s] = ismember(su_keys(:,1), validnt(ntsort));
[r, ri] = sort(s);
sortsukeys = su_keys(ri,:);
supspikes = cellfun(@(x,y) [...
    repmat(x(1),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
    repmat(x(2),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
    spikes{day}{epoch}{x(1)}{x(2)}.data ...
    repmat(y,length(spikes{day}{epoch}{x(1)}{x(2)}.data),1)], ...
    num2cell(sortsukeys,2), ...
    num2cell(1:size(sortsukeys,1))', 'un', 0);
sustack = cell2mat(supspikes(find(cell2mat(cellfun(@(x) ~isempty(x),supspikes,'un',0)))));
sustack = sortrows(sustack,3);

% Process MU into mustack (ntrode, cluster, spiketime)
mu_keys = evaluatefilter(cellinfo{day}{epoch}, 'isequal($tags,{''mua''})');
mupspikes = cellfun(@(x) [...
    repmat(x(1),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
    repmat(x(2),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
    spikes{day}{epoch}{x(1)}{x(2)}.data ...
    repmat(find(validnt(ntsort)==x(1)),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1)], ...
    num2cell(mu_keys,2), 'un', 0);
mustack = cell2mat(mupspikes(find(cell2mat(cellfun(@(x) ~isempty(x)&&size(x,2)==4, ...
    mupspikes, 'un', 0)))));
if ~isempty(mustack)
    mustack = sortrows(mustack,3);
end

%% proces to get lick bout intervals
lickboutvec = getLickBout(andef{2}, animal, [day epoch]);
bouts = evaluatefilter2(lickboutvec, '($lickBout == 1)');
nobouts = evaluatefilter2(lickboutvec, '(($lickBout == 0) & ($velocity < 1) & ($timeFromLick > .5))');
boutIntvs = vec2list(bouts{day}{epoch}(:,2), bouts{day}{epoch}(:,1));
noboutIntvs = vec2list(nobouts{day}{epoch}(:,2), nobouts{day}{epoch}(:,1));

%% plot
Pp=load_plotting_params({'defaults','dfa_plotDataChunks'});
for ts = 1:length(timeSplits(:,1))
    tstart = timeSplits(ts,1);
    tend = timeSplits(ts,2);
    sprtit = sprintf('DataChunk %s %d %d %1.f-%1.f', animal,day,epoch,tstart,tend);
    strsave = strrep(sprtit,' ', '_');
    if centerEvents
        outdir = sprintf('%s/dfa_plotDataChunk_ripcentered/%s/', pconf.andef{4},animal); 
    else
        outdir = sprintf('%s/dfa_plotDataChunk/%s/', pconf.andef{4},animal); end
    if skipExist && exist([outdir strsave], 'file'); fprintf('skipping, exists\n');
        continue; end
    if savefigs && ~pausefigs; close all;
        ifig = figure('Visible','off','units','normalized','position', Pp.position);
    else; ifig = figure('units','normalized','position',Pp.position); end
    set(gcf,'color','white')
    
    %% plot rip consensus power trace
    sf1 = subaxis(Pp.nrows,1,1,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
        'MarginBottom', Pp.MgBm); set(gca, 'Tag', 'ripkons');
    h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p=pan; p.Motion='horizontal';
    time = events{day}{epoch}{1}.eegtimesvec_ref';
    pwrIdx = knnsearch(time,[tstart; tend]);
    rippwrinWin = zrippwr(pwrIdx(1):pwrIdx(2));
    %         noiseinWin = zNoiseRipRatio(pwrIdx(1):pwrIdx(2));
    %         plot(time(pwrIdx(1):pwrIdx(2))',noiseinWin,'LineWidth',2,'Color',[.8 .3 .3 .5]);
    plot(time(pwrIdx(1):pwrIdx(2))', rippwrinWin,'LineWidth',1,'Color',[.3 .3 .3]);
    hold on;
    set(gca, 'YGrid', 'off', 'XGrid', 'on','TickDir','out', 'TickLength', [0.001 0]);
    ylabel('rippwr'); set(gca,'Xticklabel',[], 'Tag', 'zrippwr'); %xticks([]);
    axis tight; yl = ylim; xlim([tstart tend]); xl = xlim;
    patch([xl(1) xl(2) xl(2) xl(1)]', [yl(1) yl(1) minstdthresh minstdthresh]',...
        'k', 'linestyle', ':', 'edgecolor', 'none', 'facealpha', Pp.threshFAlpha);
    bpinWin = bpzpwr(pwrIdx(1):pwrIdx(2));
    plot(time(pwrIdx(1):pwrIdx(2))',bpinWin,'LineWidth',2,'Color',Pp.bpzrippwr);
    xticks(ceil(xl(1)):2:floor(xl(2))); hold off
    
    
    %% Plot LFP
    sf2 = subaxis(Pp.nrows,1,2:4); set(gca, 'Tag', 'lfp');
    lfpIdx = knnsearch(lfptime,[tstart; tend]);
    lfpCh = cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2))', ...
        eeg{day}{epoch}(validnt),'un',0)');
    gridmat = Yoffset*repmat((1:length(lfpCh(:,1)))', 1,length(lfpCh(1,:)));
    eeggrid = (lfpCh(ntsort,:) + gridmat)';
    plfp = plot(lfptime(lfpIdx(1):lfpIdx(2))', eeggrid,'LineWidth',1);
    ylabel('LFP'); hold on; axis tight; yl = ylim; xlim([tstart tend]); xl = xlim;
    set(gca,'TickDir','out'); set(plfp, {'color'}, icolors);
    set(gca, 'YGrid', 'off', 'XGrid', 'on','Xticklabel',[],'TickLength', [0.001 0], ...
        'Tag', 'ndata'); %
    ntrodeOffsets = gridmat(:,1);
    yticklabels(validnt(ntsort)); yticks([ntrodeOffsets]);
    set(gca,'fontsize',6);
    xticks(ceil(xl(1)):2:floor(xl(2)))
    
    %% Plot ripple power lfp
    sf3 = subaxis(Pp.nrows,1,5:6);
    lfpCh = zscore(double(cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2),3),...
        ripple{day}{epoch}(validnt),'un',0))))';
    ntrodeOffsets = [length(lfpCh(:,1)):-1:1]-.5;
    im = imagesc(lfptime(lfpIdx(1):lfpIdx(2))',ntrodeOffsets,lfpCh(ntsort,:));
    hold on; colormap(sf3,Pp.rippwrcmap); ylabel('ripPWR');
    axis tight; yl = ylim; xlim([tstart tend]); xl = xlim; set(gca,'Xticklabel',[]); %xticks([]);
    yticks(fliplr(ntrodeOffsets)); yticklabels(fliplr(validnt(ntsort)));
    set(gca,'TickDir','out','fontsize',4, 'YGrid', 'off', 'XGrid', 'off', ...
        'TickLength', [0.001 0], 'Tag', 'image');
    ai = find(diff(areaids));
    line(xl, [ai'; ai'], 'Color',[1 1 1 Pp.areaSepAlpha],'LineWidth', Pp.areaSepWidth);
    line([xl(1)-.01 xl(1)-.01], [yl(2) ai(2)], 'color',Pp.smclr, 'linewidth', 4);
    line([xl(1)-.01 xl(1)-.01], [ai(2) ai(1)], 'color',Pp.dmclr, 'linewidth', 4); hold off;
    
    %% Plot theta power lfp
    sf4= subaxis(Pp.nrows,1,7:8);
    lfpCh = zscore(double(cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2),3),...
        theta{day}{epoch}(validnt),'un',0))))';
    ntrodeOffsets = [length(lfpCh(:,1)):-1:1]-.5;
    im = imagesc(lfptime(lfpIdx(1):lfpIdx(2))',ntrodeOffsets,lfpCh(ntsort,:));
    hold on; colormap(sf4,Pp.thetapwrcmap); ylabel('thetaPWR');
    axis tight; yl = ylim; xlim([tstart tend]); xl = xlim; set(gca,'Xticklabel',[]); %xticks([]);
    yticks(fliplr(ntrodeOffsets)); yticklabels(fliplr(validnt(ntsort)));
    set(gca,'fontsize',4, 'TickDir','out', 'YGrid', 'off', 'XGrid', 'off', ...
        'TickLength', [0.001 0], 'Tag', 'image');
    line(xl, [ai'; ai'], 'Color',[1 1 1 Pp.areaSepAlpha],'LineWidth', Pp.areaSepWidth);
    line([xl(1)-.01 xl(1)-.01],[yl(2) ai(2)],'color',Pp.smclr,'linewidth', 4);
    line([xl(1)-.01 xl(1)-.01],[ai(2) ai(1)],'color',Pp.dmclr,'linewidth', 4); hold off;
    
    %% Plot theta phase
    sf5= subaxis(Pp.nrows,1,9:10);
    lfpCh = double(cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2),2),...
        theta{day}{epoch}(validnt),'un',0)))';
    ntrodeOffsets = [length(lfpCh(:,1)):-1:1]-.5;
    im = imagesc(lfptime(lfpIdx(1):lfpIdx(2))',ntrodeOffsets,lfpCh(ntsort,:));
    imagesc(lfptime(lfpIdx(1):lfpIdx(2))',ntrodeOffsets,flipud(lfpCh(ntsort,:)));
    hold on; colormap(sf5, Pp.phaseCmap); ylabel('thetaPHASE');
    axis tight; yl = ylim; xlim([tstart tend]); xl = xlim; set(gca,'Xticklabel',[]); %xticks([]);
    yticks(fliplr(ntrodeOffsets)); yticklabels(fliplr(validnt(ntsort)));
    set(gca,'fontsize',4,'TickDir','out', 'YGrid', 'off', 'XGrid', 'off', ...
        'TickLength', [0.001 0], 'Tag', 'image');
    line(xl, [ai'; ai'], 'Color',[1 1 1 Pp.areaSepAlpha],'LineWidth', Pp.areaSepWidth);
    line([xl(1)-.01 xl(1)-.01], [yl(2) ai(2)],'color',Pp.smclr,'linewidth', 4);
    line([xl(1)-.01 xl(1)-.01], [ai(2) ai(1)],'color',Pp.dmclr,'linewidth', 4);hold off;
    
    %% Plot SU spikes
    suCh = sustack(all([sustack(:,3)>tstart sustack(:,3)<tend],2),:);
    if ~isempty(suCh)
        sf6= subaxis(Pp.nrows,1,11:12);
        uclr = zeros(size(suCh,1),3);
        uclr(ismember(suCh(:,1),c1Nt),:)=repmat(Pp.c1clr,sum(ismember(suCh(:,1),c1Nt)),1);
        uclr(ismember(suCh(:,1),dmNt),:)=repmat(Pp.dmclr,sum(ismember(suCh(:,1),dmNt)),1);
        uclr(ismember(suCh(:,1),smNt),:)=repmat(Pp.smclr,sum(ismember(suCh(:,1),smNt)),1);
        dx = scatter(suCh(:,3), suCh(:,4),Pp.SUsize,uclr,'filled'); dx.Marker = 'd';
        ylabel('SU'); set(gca, 'TickDir','out', 'YGrid', 'off', 'XGrid', 'on', ...
            'fontsize', 4,'TickLength', [0.001 0], 'Tag', 'ndata', 'Xticklabel',[]);
        axis tight; yl = ylim; ylim([yl(1)-2 yl(2)+2]);
        xlim([tstart tend]); xl = xlim; xticks(ceil(xl(1)):2:floor(xl(2)));
        [yt, ydx] = unique(suCh(:,4)); yticks(yt); yticklabels(suCh(ydx,1));
    end
    
    %% Plot MU spikes
    if any(mustack)
        muCh = mustack(all([mustack(:,3)>tstart mustack(:,3)<tend],2),:);
    if ~isempty(muCh)
        sf7= subaxis(Pp.nrows,1,13);
        uclr = zeros(size(muCh,1),3);
        uclr(ismember(muCh(:,1),c1Nt),:)=repmat(Pp.c1clr,sum(ismember(muCh(:,1),c1Nt)),1);
        uclr(ismember(muCh(:,1),dmNt),:)=repmat(Pp.dmclr,sum(ismember(muCh(:,1),dmNt)),1);
        uclr(ismember(muCh(:,1),smNt),:)=repmat(Pp.smclr,sum(ismember(muCh(:,1),smNt)),1);
        dx = scatter(muCh(:,3), muCh(:,4),1,uclr,'filled'); ylabel('MU');
        ylabel('SU'); set(gca,'TickDir','out','YGrid','off','XGrid','on','fontsize', 4,...
            'TickLength', [0.001 0],'Tag', 'ndata', 'Xticklabel',[]);
        axis tight; yl = ylim; ylim([yl(1)-2 yl(2)+2]); %
        xlim([tstart tend]); xl = xlim; xticks(ceil(xl(1)):1:floor(xl(2)))
        [yt, ydx] = unique(muCh(:,4)); yticks(yt); yticklabels(muCh(ydx,1));
    end; end

    %% plot velocity, pos, hd
    sf8= subaxis(Pp.nrows,1,14);
    velcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'vel-loess'));
    xcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'x-loess'));
    ycol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'y-loess'));
    dircol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'dir-loess'));
    posIdx = knnsearch(postime(:,1),[tstart; tend]);
    vel = pos{day}{epoch}.data(posIdx(1):posIdx(2),velcol);
    veltime = postime(posIdx(1):posIdx(2));
    plot(veltime, vel); hold on; axis tight; % just to get ylims
    ylabel('Vel2D'); axis tight; set(gca,'YGrid','off','XGrid','on','TickDir','out','TickLength', [0.001 0]);
    yl = ylim; xlim([tstart tend]); xl = xlim; xticks(ceil(xl(1)):2:floor(xl(2)))
    velinds = vel < maxvelocity; v = vel(velinds); v(v<4) = 1; v(v>=4) = 0; set(gca,'Xticklabel',[]); %
    stem(veltime(velinds),v*yl(2),'k','filled','Marker', 'none', 'linewidth',2,'color',[.8 .8 .8 .2]);
    plot(veltime,vel,'LineWidth',2,'Color', 'k', 'linestyle', ':'); hold off
    
    %% lin pos
    sf9= subaxis(Pp.nrows,1,15);
    lpIdx = knnsearch(linpostime(:,1),[tstart; tend]);
    lpClr = linpos{day}{epoch}.statematrix.segmentIndex(lpIdx(1):lpIdx(2));
    linpostimeCh = linpostime(lpIdx(1):lpIdx(2));
    lindistCh = lindist1D(lpIdx(1):lpIdx(2));
    scatter(linpostimeCh,lindistCh,7,lpClr,'filled', 'marker', 'd'); hold on;
    colormap(Pp.linposcmap); ylabel('linPos'); xlabel('Time s');
    axis tight; xlim([tstart tend]); xl = xlim; xticks(ceil(xl(1)):2:floor(xl(2)))
    set(gca,'TickDir','out','YGrid','off', 'XGrid', 'on','TickLength', [0.001 0], ...
        'Tag', 'linpos'); %
    for iseg = 1:length(segSE(:,1))
        patch(sf9, [xl(1) xl(2) xl(2) xl(1)], ...
            [segSE(iseg,1) segSE(iseg,1) segSE(iseg,2) segSE(iseg,2)], iseg,...
            'FaceAlpha',Pp.linposSegAlpha,'edgecolor','none'); end; hold off;
    
    %% plot dio stems
    isinput = cellfun(@(x) isequal(x.input,1), DIO{day}{epoch}, 'un', 1);
    dioID = cellfun(@(x) str2double(regexp(x.original_id,'\d*','Match')), ...
        DIO{day}{epoch}, 'un', 1);
    outputdios = task{day}{epoch}.outputdio;
    inputdios = task{day}{epoch}.inputdio;
    outdioIdx = find(all([ismember(dioID, outputdios)' ~isinput'],2));
    indioIdx = find(all([ismember(dioID, inputdios)' isinput'],2));
    
    IOdioCh =[indioIdx; outdioIdx];
    dDF = []; diosDFfields = {'chan', 'startTime', 'endTime', 'isOutputWell'};
    for c = 1:length(IOdioCh) % collect input dios
        ch = IOdioCh(c); diotimes = [];
        % valid events are bool flips..
        dioChIdx = all([DIO{day}{epoch}{ch}.times>tstart ...
            DIO{day}{epoch}{ch}.times<tend],2);
        if sum(dioChIdx) == 0; continue; end
        dioTimeCh = double(DIO{day}{epoch}{ch}.times(dioChIdx));
        dioValsCh = double(DIO{day}{epoch}{ch}.values(dioChIdx));
        % dios are Up (1) to Down (0)
        while ~isempty(dioValsCh) && dioValsCh(1) == 0
            dioValsCh(1) = []; dioTimeCh(1) = []; end
        while ~isempty(dioValsCh) && dioValsCh(end) == 1
            dioValsCh(end) = []; dioTimeCh(end) = []; end
        if ~isempty(dioValsCh)
            dioVd = [1; find(abs(diff(dioValsCh)))+1];
            diotimes = dioTimeCh(dioVd);
            dioStEnd = [diotimes(1:2:end) diotimes(2:2:end)];
            isout = any(ismember(outputdios, ch));
            dDF = [dDF; ch*ones(length(dioStEnd(:,1)),1) dioStEnd ...
                isout*ones(length(dioStEnd(:,1)),1)];
        else; continue; end; end
    if ~isempty(dDF)
        c = colorPicker(cellstr(string(dDF(:,1)))); % color by dio chan id
        ax = get(ifig, 'children');
        for ix = 1:length(ax)
            if strcmp(ax(ix).Tag, 'image')
                continue; end
            try
%                 disp(ax(ix).Tag)
                if strcmp(ax(ix).Tag, 'zrippwr')
                    yl = ylim(ax(ix)); hold(ax(ix),'on')
                    opts = dDF(:,1)>3;
                    scatter(ax(ix), dDF(find(opts),2),...
                        repmat(yl(2),size(dDF(find(opts),1),1),1),50,...
                        dDF(find(opts),1)>3, 'o', 'filled');
%                     scatter(ax(ix), dDF(find(~opts),2), ...
%                         repmat(yl(2),size(dDF(find(~opts),1),1),1),30,...
%                         dDF(find(~opts),1)>3, '+');
                    hold(ax(ix),'off')
                end
                l = line(ax(ix),dDF(:,[2 2])',ylim(ax(ix)), 'linewidth',.5);
                if length(c(:,1)) > 1
                    set(l,{'color'},c);
                else
                    set(l,'color',c);
                end
                if strcmp(ax(ix).Tag, 'ndata')
                    set(ax(ix),'children',flipud(get(ax(ix),'children')))
                end
            catch
                continue
            end; end; end
    %% plot ripple patches
    swrInWinIdx = all([eventTime(:,1)>tstart eventTime(:,2)<tend],2);
    if ~any(swrInWinIdx); fprintf('no rips in win %1.f : %1.f\n', tstart, tend);
        if skipNoRipChunks; fprintf('skipping\n');
            continue;
        end
    end
    if any(swrInWinIdx)
        swrInWin = eventTime(swrInWinIdx,:);
        xs = swrInWin(:,1);
        xe = swrInWin(:,2);
        ax = get(ifig, 'children');
        for ix = 1:length(ax)
            if strcmp(ax(ix).Tag, 'image')
                continue;
            end
            try
                %                 disp(ax(ix).Tag)
                yl = ylim(ax(ix));
                if strcmp(ax(ix).Tag, 'zrippwr')
                    hold(ax(ix),'on')
                    s = scatter(ax(ix), xs,repmat(yl(2),size(xs,1),1),30,'filled',...
                        'marker', 'd', 'markerfacecolor', 'y', 'markeredgecolor', 'k');
                    hold(ax(ix),'off')
                end
                patch(ax(ix), [xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
                    length(xe)),Pp.ripclr, 'FaceAlpha',Pp.ripFAlpha, 'edgecolor','none');
            catch
                continue
            end; end; end
    
    %% plot lick bout patches
    evs.bouts = boutIntvs;
    evs.nobouts = noboutIntvs;
    fns = fieldnames(evs);
    for ev = 1:length(fns)
        fieldn = fns{ev};
        eventIntervals = evs.(fieldn);
        eventsInWinIdx = all([eventIntervals(:,1)>tstart eventIntervals(:,2)<tend],2);
        if any(eventsInWinIdx)
            eventInWin = eventIntervals(eventsInWinIdx,:);
            xs = eventInWin(:,1);
            xe = eventInWin(:,2);
            ax = get(ifig, 'children');
            for ix = 1:length(ax)
                if strcmp(ax(ix).Tag, 'image')
                    continue;
                end
                try
                    hold(ax(ix),'on')
                    yl = ylim(ax(ix));
                    if strcmp(ax(ix).Tag, 'zrippwr')
                        s = scatter(ax(ix), xs,repmat(yl(2),size(xs,1),1),30,'filled',...
                            'marker','*', 'markerfacecolor', Pp.([fieldn 'clr']), 'markeredgecolor','k');
                    end
                    patch(ax(ix), [xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
                        length(xe)),Pp.([fieldn 'clr']), 'FaceAlpha',Pp.([fieldn 'alpha']), 'edgecolor','none');
                    hold(ax(ix),'off')
                catch
                    continue
                end; end
        else
            fprintf('no %s in win %1.f : %1.f\n', fns{ev}, tstart, tend);
            if skipNoEventChunks; fprintf('skipping\n');
                continue;
            end
        end
    end
    %% ----- link x axis -----
    allAxesInFigure = findall(ifig,'type','axes'); linkaxes(allAxesInFigure, 'x');
    % ---- super axis -----
    sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
    iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
    set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment','center','FontSize',12);
    h = get(gcf,'Children');
    set(gcf,'Children',flip(h)); % put super axis at bottom of axis stack. allows for zoom
    % ---- pause, save figs ----
    if pausefigs; pause; end
    if savefigs; save_figure(outdir, sprtit, 'savefigas', savefigas); end
end
end