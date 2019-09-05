

%{

Change this toeither plot all datachunks centered on each ripple.. 
OR plot all the 10s timechunks, regaredless of whether there were ripples
during it or not
%}

create_filter =1;
run_ff = 0;
save_ffdata = 0;
load_ffdata = 0;

loadEpoch = 1;
processData = 1;

plotfigs = 1;
skipNoRipChunks = 1;
skipExist = 1;
centerEvents = 1;
centerOn = 'ca1ripplekons';
savefigs = 0;
pausefigs = 1;
%%
pconf = paramconfig;
Fp.animals = {'D10'};%, 'JZ1', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
Fp.filtfunction = 'dfa_plotDataChunks';
Fp.add_params = {'wtrackdays','excludeNoise','excludePriorFirstWell','<4cm/s', ...
    'excludeAfterLastWell', 'valid_ntrodes'};
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
day = 5;
epoch = 4;
splitSize = 10; % seconds
Yoffset = 500;
%% FF
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
       'excludetime', Fp.timefilter,'iterator',Fp.iterator,'eegtetrodes',Fp.tetfilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff; F = runfilter(F); for d = 1:length(F); F(d).datafilter_params = Fp; end; end
if save_ffdata
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s', Fp.epochEnvironment)); end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s', Fp.epochEnvironment)); end
%% run one epoch for dev
if loadEpoch
        % replaces the iterator to load a single epoch
    animal = Fp.animals{1};
    animdef = animaldef(animal);
    events = loaddatastruct(animdef{2}, animal, 'ca1rippleskons', day);
    noisekons = loaddatastruct(animdef{2}, animal, 'ca1noisekons', day);
    eeg = loadeegstruct(animdef{2}, animal, 'eeg', day, epoch, F.eegdata{1}{2});
    ripple = loadeegstruct(animdef{2}, animal, 'ripple', day, epoch, F.eegdata{1}{2});
    theta = loadeegstruct(animdef{2}, animal, 'theta', day, epoch, F.eegdata{1}{2});
    pos = loaddatastruct(animdef{2}, animal, 'pos', day);
    linpos = loaddatastruct(animdef{2}, animal, 'linpos', day);
    spikes = loaddatastruct(animdef{2}, animal, 'spikes', day);
    cellinfo = loaddatastruct(animdef{2}, animal, 'cellinfo', day);
    tetinfo = loaddatastruct(animdef{2}, animal, 'tetinfo', day);
    dio = loaddatastruct(animdef{2}, animal, 'DIO', day);
    task = loaddatastruct(animdef{2}, animal, 'task', day);
end
if processData

    %% this is where the dfa starts
%     noiseEvents=load_data('filterframework','noiseEvents', Fp.animals, 'animpos', 0);
    try
        excludeIntervals = F.excludetime{1}{ismember(F.epochs{1},[day epoch], 'rows')};
    catch
        error('no data for %d %d\n', day, epoch);
    end
    eventTime = [events{day}{epoch}{1}.starttime events{day}{epoch}{1}.endtime];
    eventTime = eventTime(~isExcluded(eventTime(:,1),excludeIntervals),:);
    % data intervals
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
        timeSplits = [startTime:splitSize:endTime]';
        timeSplits = [timeSplits(1:end-1) timeSplits(2:end)]; % start end
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
    lindist1D(outerIdx) = lindist1D(outerIdx)+innermax;
    segStartEnd = []; c = 0;
    % get bounds of lindist segments
    for iseg = 1:length(unique(linpos{day}{epoch}.statematrix.segmentIndex))
        segstart = min(lindist1D(linpos{day}{epoch}.statematrix.segmentIndex == iseg));
        segend = max(lindist1D(linpos{day}{epoch}.statematrix.segmentIndex == iseg));
        if ~isempty(segstart)
            c = c+1;
            segStartEnd(c,1) = segstart;
            segStartEnd(c,2) = segend;
        end
    end
    % rip zscpew
    zrippwr = ((events{day}{epoch}{1}.powertrace-nanmean(events{day}{epoch}{1}.powertrace))...
        ./nanstd(events{day}{epoch}{1}.powertrace))';
    zNoiseRipRatio = ((noisekons{day}{epoch}{1}.powertrace-nanmean(noisekons{day}{epoch}{1}.powertrace))...
        ./nanstd(noisekons{day}{epoch}{1}.powertrace))';
    % get nt's for each area
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
        'mec', 'supf'; 'mec', 'nsupf'};
%         'ref', 'ca1'; 'ref', 'mec'};
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
    icolors = colorPicker(arTag, sbTag); icolors = icolors(ntsort);
    
    % GET SU (ntrode, cluster, spiketime, cellnum) 
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
    % Get MU (ntrode, cluster, spiketime) 
    mu_keys = evaluatefilter(cellinfo{day}{epoch}, 'isequal($tags,{''mua''})');
    mupspikes = cellfun(@(x) [...
        repmat(x(1),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
        repmat(x(2),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
        spikes{day}{epoch}{x(1)}{x(2)}.data ...
        repmat(find(validnt(ntsort)==x(1)),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1)], ...
        num2cell(mu_keys,2), 'un', 0);
    mustack = cell2mat(mupspikes(find(cell2mat(cellfun(@(x) ~isempty(x)&&size(x,2)==4, ...
        mupspikes, 'un', 0)))));
    mustack = sortrows(mustack,3);
end
%% plot
if plotfigs
    Pp=load_plotting_params({'defaults','dataExplore'});
    nrows = Pp.nrows;
    for ts = 1:length(timeSplits(:,1))
        tstart = timeSplits(ts,1);
        tend = timeSplits(ts,2);
        sprtit = sprintf('DataChunk %s %d %d %1.f-%1.f', animal, day, epoch, tstart, tend);
        strsave = strrep(sprtit,' ', '_');
        outdir = sprintf('%s/dfa_plotDataChunk/%s/',pconf.andef{4},animal);
        if skipExist && exist([outdir strsave], 'file')
            continue
        end
        % get ripple times in window
        swrInWinIdx = all([eventTime(:,1)>tstart eventTime(:,2)<tend],2);
        if ~any(swrInWinIdx)
            fprintf(sprintf('no rips in win %1.f : %1.f\n', tstart, tend)); 
            if skipNoRipChunks
                continue;
            end
        end
        if savefigs && ~pausefigs; close all;
            ifig = figure('Visible','off','units','normalized','position', Pp.position);
        else; ifig = figure('units','normalized','position',Pp.position); end
        set(gcf,'color','white')
        %% plot rip consensus power trace
        sf1 = subaxis(nrows,1,1,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
            Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
        time = events{day}{epoch}{1}.eegtimesvec_ref';
        pwrIdx = knnsearch(time,[tstart; tend]);
        rippwrinWin = zrippwr(pwrIdx(1):pwrIdx(2));
        noiseinWin = zNoiseRipRatio(pwrIdx(1):pwrIdx(2));
        plot(time(pwrIdx(1):pwrIdx(2))',noiseinWin,'LineWidth',2,'Color',[.8 .3 .3 .5]); hold on; 
        plot(time(pwrIdx(1):pwrIdx(2))', rippwrinWin,'LineWidth',2,'Color',[.3 .3 .3]);
        h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p = pan; p.Motion = 'horizontal';
        set(gca, 'YGrid', 'off', 'XGrid', 'on'); ylabel('rippwr'); xticks([]); 
        set(gca,'TickDir','out'); axis tight; yl = ylim; xl = xlim; 
        patch([xl(1) xl(2) xl(2) xl(1)]', [yl(1) yl(1) Fp.minstdthresh Fp.minstdthresh]',...
            'k', 'linestyle', ':', 'edgecolor', 'none', 'facealpha', Pp.threshFAlpha);

        %% Plot LFP
        sf2 = subaxis(nrows,1,2:4);%, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
%             Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
        lfpIdx = knnsearch(lfptime,[tstart; tend]);
        lfpCh = cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2))', ...
            eeg{day}{epoch}(validnt),'un',0)');
        gridmat = Yoffset*repmat((1:length(lfpCh(:,1)))', 1,length(lfpCh(1,:)));
        eeggrid = (lfpCh(ntsort,:) + gridmat)';
        plfp = plot(lfptime(lfpIdx(1):lfpIdx(2))', eeggrid,'LineWidth',1);
        h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p = pan; p.Motion = 'horizontal'; 
        ylabel([]); hold on; axis tight; yl = ylim; xticks([]);
        set(gca,'TickDir','out'); set(plfp, {'color'}, icolors); 
        set(gca, 'YGrid', 'off', 'XGrid', 'on');
        a_sA = cellfun(@(x,y,z) sprintf('%s %s %d',x,y,z),arTag(ntsort), ...
            sbTag(ntsort),num2cell(validnt(ntsort)'),'un', 0);
        ntrodeOffsets = gridmat(:,1);
        yticklabels(a_sA); yticks([ntrodeOffsets]); ytickangle(0);
        a = get(gca,'YTickLabel'); set(gca,'YTickLabel',a,'fontsize',6); 
        

        %% Plot ripple power lfp
        sf3 = subaxis(nrows,1,5:6);%, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
%             'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
%             'MarginBottom', Pp.MgBm);
        % gath envelope magnitude of riplfp from this window and Vert Offset for plotting
        lfpCh = double(cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2),3),...
            ripple{day}{epoch}(validnt),'un',0)))';
        im = imagesc(lfptime(lfpIdx(1):lfpIdx(2))', [1:length(lfpCh(:,1))]-.5, ...
            flipud(lfpCh(ntsort,:)));
        colormap(sf3,Pp.rippwrcmap); 
        h=zoom;h.Motion='horizontal';h.Enable='on';p=pan;p.Motion='horizontal'; 
        ylabel('ripLFP'); axis tight; yl = ylim; xticks([]); set(gca,'TickDir','out');
        yticks(ntrodeOffsets); yticklabels(a_sA); ai = find(diff(areaids)); hold on;
        line(xl, [ai'; ai'], 'Color',[1 1 1 Pp.areaSepAlpha],'LineWidth', Pp.areaSepWidth);
        line([xl(1)-.01 xl(1)-.01], [yl(2) ai(2)], 'color',Pp.smclr, 'linewidth', 4);
        line([xl(1)-.01 xl(1)-.01], [ai(2) ai(1)], 'color',Pp.dmclr, 'linewidth', 4); hold off;
        
       %% Plot theta power lfp
        sf4= subaxis(nrows,1,7:8);% 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
%             'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
%             'MarginBottom', Pp.MgBm);
        % envelope magnitude of thetalfp in window with vertical offset
        lfpCh = double(cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2),3),...
            theta{day}{epoch}(validnt),'un',0)))';
        im = imagesc(lfptime(lfpIdx(1):lfpIdx(2))', [1:length(lfpCh(:,1))]-.5,...
            flipud(lfpCh(ntsort,:)));
        colormap(sf4,Pp.thetapwrcmap);
        h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p = pan; p.Motion = 'horizontal'; 
        ylabel('thetaLFP'); hold on; axis tight; yl = ylim; xticks([]); set(gca,'TickDir','out');
        yticks([ntrodeOffsets]); yticklabels(a_sA);
        line(xl, [ai'; ai'], 'Color',[1 1 1 Pp.areaSepAlpha],'LineWidth', Pp.areaSepWidth);
        line([xl(1)-.01 xl(1)-.01], [yl(2) ai(2)],'color',Pp.smclr,'linewidth', 4);
        line([xl(1)-.01 xl(1)-.01], [ai(2) ai(1)],'color',Pp.dmclr,'linewidth', 4);hold off;
        
        %% Plot theta phase
        sf5= subaxis(nrows,1,9:10); %, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
%             'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
%             'MarginBottom', Pp.MgBm);
        lfpCh = double(cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2),2),...
            theta{day}{epoch}(validnt),'un',0)))';
        imagesc(lfptime(lfpIdx(1):lfpIdx(2))',1:length(lfpCh(:,1))-.5,...
            flipud(lfpCh(ntsort,:)));
        colormap(sf5, Pp.phaseCmap);
        h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p = pan; p.Motion = 'horizontal'; 
        ylabel('thetaLFP'); axis tight; yl = ylim; xticks([]); set(gca,'TickDir','out'); 
        yticks([ntrodeOffsets]); yticklabels(a_sA); hold on;
        line(xl, [ai'; ai'], 'Color',[1 1 1 Pp.areaSepAlpha],'LineWidth', Pp.areaSepWidth);
        line([xl(1)-.01 xl(1)-.01], [yl(2) ai(2)],'color',Pp.smclr, 'linewidth', 4);
        line([xl(1)-.01 xl(1)-.01], [ai(2) ai(1)],'color',Pp.dmclr,'linewidth',4);hold off;
        
        %% Plot SU spikes
        sf6= subaxis(nrows,1,11:12);%,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
%             Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
        suCh = sustack(all([sustack(:,3)>tstart sustack(:,3)<tend],2),:);
        uclr = zeros(size(suCh,1),3);
        uclr(ismember(suCh(:,1),c1Nt),:)=repmat(Pp.c1clr,sum(ismember(suCh(:,1),c1Nt)),1);
        uclr(ismember(suCh(:,1),dmNt),:)=repmat(Pp.dmclr,sum(ismember(suCh(:,1),dmNt)),1);
        uclr(ismember(suCh(:,1),smNt),:)=repmat(Pp.smclr,sum(ismember(suCh(:,1),smNt)),1);
        a = scatter(suCh(:,3), suCh(:,4),Pp.SUsize,uclr,'filled');
        a.Marker = 'd'; set(gca,'TickDir','out'); set(gca, 'YGrid', 'off', 'XGrid', 'on');
        ylabel('suSpikes'); axis tight; yl = ylim; ylim([yl(1)-2 yl(2)+2]); xticks([]); 

        %% Plot MU spikes
        sf7= subaxis(nrows,1,13);%,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
%             Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
        muCh = mustack(all([mustack(:,3)>tstart mustack(:,3)<tend],2),:);
        uclr = zeros(size(muCh,1),3);
        uclr(ismember(muCh(:,1),c1Nt),:)=repmat(Pp.c1clr,sum(ismember(muCh(:,1),c1Nt)),1);
        uclr(ismember(muCh(:,1),dmNt),:)=repmat(Pp.dmclr,sum(ismember(muCh(:,1),dmNt)),1);
        uclr(ismember(muCh(:,1),smNt),:)=repmat(Pp.smclr,sum(ismember(muCh(:,1),smNt)),1);
        a = scatter(muCh(:,3), muCh(:,4),1,uclr,'filled');
        ylabel('muSpikes'); axis tight; yl = ylim; ylim([yl(1)-2 yl(2)+2]); xticks([]); 
        set(gca,'TickDir','out'); set(gca, 'YGrid', 'off', 'XGrid', 'on');
        
        %% plot velocity, pos, hd
        sf8= subaxis(nrows,1,14);%,'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
%             Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
        velcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'vel-loess'));
        xcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'x-loess'));
        ycol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'y-loess'));
        dircol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'dir-loess'));
        posIdx = knnsearch(postime(:,1),[tstart; tend]);
        vel = pos{day}{epoch}.data(posIdx(1):posIdx(2),velcol);
        veltime = postime(posIdx(1):posIdx(2));
        plot(veltime, vel, 'LineWidth', 2, 'Color', 'k', 'linestyle', '-');
        ylabel('2D Vel'); axis tight; xticks([]); set(gca,'YGrid','off','XGrid','on');
        set(gca,'TickDir','out'); yl = ylim; xl = xlim;
        velinds = vel < Fp.maxvelocity; v = vel(velinds); v(v<4) = 1; v(v>=4) = 0; hold on;
        stem(veltime(velinds),v*yl(2),'.k','filled','linewidth',5,'color',[.8 .8 .8 .2]);
        plot(veltime,vel,'LineWidth',2,'Color', 'k', 'linestyle', '-'); hold off
        
        %% lin pos
        sf9= subaxis(nrows,1,15, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, 'MarginLeft', ...
            Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);

        linposIdx = knnsearch(linpostime(:,1),[tstart; tend]);
        linpostimeCh = linpostime(linposIdx(1):linposIdx(2));
        lindistCh = lindist1D(linposIdx(1):linposIdx(2));
        plot(linpostimeCh,lindistCh,'linewidth',3,'color',[.8 .8 .8]); %background pos
        xlabel('Time s'); ylabel('Lindist'); hold on; 
        set(gca,'TickDir','out'); set(gca, 'YGrid', 'on', 'XGrid', 'on');
        axis tight; yl = ylim;
%         ylim([0 200]); 
%         segclr = lines(5);
%         clr1 = [.6 .5 .2; .8 .5 .8; .6 0 .6; 0 0.4980 0.3255; 0 0.8 0.5216];
%         segmentnum = linpos{day}{epoch}.statematrix.segmentIndex(linposIdx(1):linposIdx(2));
        for iseg = 1:length(segStartEnd(:,1))
            patch(sf9, [xl(1) xl(1)+.1 xl(1)+.1 xl(1)], ...
                    [segStartEnd(iseg,1) segStartEnd(iseg,1) segStartEnd(iseg,2) segStartEnd(iseg,2)], iseg,'FaceAlpha',.9,'edgecolor','none');
        end
%         line([xl(1) xl(1)-2], [ai(2) ai(1)],'color',Pp.dmclr,'linewidth',4);hold off;
        %% plot dios over linpos
        outputdios = task{day}{epoch}.outputdio;
        IOdioCh =[task{day}{epoch}.inputdio task{day}{epoch}.outputdio];
        diosDF = [];
        diosDFfields = {'chan', 'startTime', 'endTime', 'isOutputWell'};
        % for each input/output dio, plot patch 
        for c = 1:length(IOdioCh) % collect input dios
            ch = IOdioCh(c); diotimes = [];
            % valid events are bool flips.. guards against crosstalk
            dioChIdx = all([dio{day}{epoch}{ch}.times>tstart dio{day}{epoch}{ch}.times<tend],2);
            if sum(dioChIdx) == 0; continue; end
            dioTimeCh = double(dio{day}{epoch}{ch}.times(dioChIdx));
            dioValsCh = double(dio{day}{epoch}{ch}.values(dioChIdx));
            dioVd = find([1; abs(diff(dioValsCh))]);
            while dio{day}{epoch}{ch}.values(dioVd(1)) == 0
                dioVd(1) = []; end
            while dio{day}{epoch}{ch}.values(dioVd(end)) == 1
                dioVd(end) = []; end 
            diotimes = dioTimeCh(dioVd);
            dioStEnd = [diotimes(1:2:end) diotimes(2:2:end)];
            isout = any(ismember(outputdios, ch));
            diosDF = [diosDF; ch*ones(length(dioStEnd(:,1)),1) dioStEnd ...
                isout*ones(length(dioStEnd(:,1)),1)];
        end
            a = [diosDF(:,2:3)'; diosDF(:,[3 2])'];
            b = repmat([yl(1) yl(1) yl(2) yl(2)]',1, length(diosDF(:,1)));
            patch(sf9,a,b,diosDF(:,1),'FaceAlpha',Pp.dioFaceAlpha, 'edgealpha', Pp.dioEdgeAlpha);
            hold off
       %% plot ripple patches
       if any(swrInWinIdx)
            swrInWin = eventTime(swrInWinIdx,:);
            xs = swrInWin(:,1);
            xe = swrInWin(:,2);
            yl = ylim(sf1);
            patch(sf1, [xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
                length(xe)),Pp.ripclr, 'FaceAlpha',Pp.ripFAlpha, 'edgecolor','none');
            yl = ylim(sf2);
            patch(sf2, [xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
                length(xe)),Pp.ripclr,'FaceAlpha',Pp.ripFAlpha, 'edgecolor','none');
            yl = ylim(sf3);
            patch(sf3, [xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1, ...
                length(xe)),Pp.ripclr, 'FaceAlpha',0, 'edgecolor','w', 'edgealpha', .5);
            yl = ylim(sf4);
            patch(sf4,[xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
                length(xe)),Pp.ripclr, 'FaceAlpha',0, 'edgecolor','k', 'edgealpha', .5); 
            yl = ylim(sf5);
            patch(sf5, [xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
                length(xe)),Pp.ripclr,'FaceAlpha',Pp.ripFAlpha, 'edgecolor','b', ...
                'edgealpha', .5);
            yl = ylim(sf6);
            patch(sf6, [xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
                length(xe)),Pp.ripclr,'FaceAlpha',Pp.ripFAlpha, 'edgecolor','none');
            yl = ylim(sf7);
            patch(sf7, [xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
                length(xe)),Pp.ripclr,'FaceAlpha',Pp.ripFAlpha, 'edgecolor','none');
            yl = ylim(sf8);
            patch(sf8, [xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
                length(xe)),Pp.ripclr,'FaceAlpha',Pp.ripFAlpha, 'edgecolor','none');
            yl = ylim(sf9);
            patch(sf9, [xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
                length(xe)),Pp.ripclr,'FaceAlpha',Pp.ripFAlpha, 'edgecolor','none');
       end
       
        %% ----- link x axis -----
        allAxesInFigure = findall(ifig,'type','axes'); linkaxes(allAxesInFigure, 'x');
        % ---- super axis -----
        sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
        iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
        set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment','center','FontSize',12);
        % ---- pause, save figs ----
        if pausefigs; pause; end
        if savefigs; save_figure(outdir, sprtit); end
    end
end

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