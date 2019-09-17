

%{
 i want to have filterframework get the data instead of plotting it
- i need to run the createfilter in order to get the appropriate
excludetimes.. much like i did with the riptrig lfp..
%}

create_filter =1;
run_ff = 1;
save_ffdata = 0;
load_ffdata = 0;
processData = 0;
plotfigs = 0;
savefigs = 0;
pausefigs = 1;
%%
pconf = paramconfig;
Fp.animals = {'D10', 'JZ1', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
Fp.filtfunction = 'dfa_plotDataChunks';
Fp.add_params = {'wtrackdays','excludeNoise','excludePriorFirstWell','<4cm/s', ...
    'excludeAfterLastWell', 'valid_ntrodes'};
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
day = 2;
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
if processData
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
    
    % rip zscpew
    zrippwr = ((events{day}{epoch}{1}.powertrace-nanmean(events{day}{epoch}{1}.powertrace))...
        ./nanstd(events{day}{epoch}{1}.powertrace))';
    zNoiseRipRatio = ((noisekons{day}{epoch}{1}.powertrace-nanmean(noisekons{day}{epoch}{1}.powertrace))...
        ./nanstd(noisekons{day}{epoch}{1}.powertrace))';
    
    % get nt's for each area
    andef = animaldef(animal);
    ntinfo = loaddatastruct(andef{2}, andef{3}, 'tetinfo');    
%     reftet_keys = evaluatefilter(ntinfo, 'isequal($area, ''ref'')');
    ca1Nts = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''ca1'')');
    mecSpNts = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''mec'') && (isequal($subarea, ''supf'') || isequal($subarea, ''nsupf''))');
    mecDpNts = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''mec'') && (isequal($subarea, ''deep'') || isequal($subarea, ''ndeep''))');

    % ---- reorder the nt's based on area,subarea ----
    areaTags = cellfun(@(x) ntinfo{day}{epoch}{x}.area,num2cell(validntrodes),'un',0)';
    subareaTags = cellfun(@(x) ntinfo{day}{epoch}{x}.subarea,num2cell(validntrodes),'un',0)';
    ordernts = {'ca1', 'd'; ...
        'mec', 'deep'; 'mec', 'ndeep'; ...
        'mec', 'supf'; 'mec', 'nsupf'};
%         'ref', 'ca1'; 'ref', 'mec'};
    areaTagsNT = {};
    subareaTagsNT = {};
%     areaTagsNT(validntrodes,1) = areaTags;
%     subareaTagsNT(validntrodes,1) = subareaTags;
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
    
    % GET SU
    su_f = 'all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''),$tags,''un'',0))))';
    su_keys = evaluatefilter(cellinfo{day}{epoch}, su_f);
    mu_keys = evaluatefilter(cellinfo{day}{epoch}, 'isequal($tags,{''mua''})');
    
    su_keys = su_keys(ismember(su_keys(:,1), validntrodes),:);
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
    sustack = cell2mat(supspikes(find(cell2mat(cellfun(@(x) ~isempty(x),supspikes,'un',0)))));
    sustack = sortrows(sustack,3);
    
    % create array (ntrode, cluster, spiketime,ntsorted)
    mupspikes = cellfun(@(x) [...
        repmat(x(1),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
        repmat(x(2),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1) ...
        spikes{day}{epoch}{x(1)}{x(2)}.data ...
        repmat(find(validntrodes(ntsort)==x(1)),length(spikes{day}{epoch}{x(1)}{x(2)}.data),1)], ...
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
        patch([xl(1) xl(2) xl(2) xl(1)]', [yl(1) yl(1) Fp.minstdthresh Fp.minstdthresh]','k', ...
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
        eegstack = double(cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2),3),...
            theta{day}{epoch}(validntrodes),'un',0)))';
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
        
        %% Plot theta phase and phase diff
        sfl = subaxis(nrows,1,9:10, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
            'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
            'MarginBottom', Pp.MgBm);
        % envelope magnitude of thetalfp in window with vertical offset
        eegstack = double(cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2),2),...
            theta{day}{epoch}(validntrodes),'un',0)))';
%         eegstack = diff([eegstack(:,1) eegstack],[],2);
        im = imagesc(lfptime(lfpIdx(1):lfpIdx(2))', [1:length(eegstack(:,1))]-.5,...
            eegstack(ntsort,:));
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
        ylabel('SU spikes'); hold on; axis tight; yl = ylim; ylim([yl(1)-2 yl(2)+2]);
        yl = ylim; set(gca, 'xticklabel', []);
        set(gca,'TickDir','out'); set(gca, 'YGrid', 'off', 'XGrid', 'on');
        patch([xs'; xe'; xe'; xs'], repmat([yl(1) yl(1) yl(2) yl(2)]', 1, length(xe)),'y', ...
            'FaceAlpha',.15, 'edgecolor','none'); hold off
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
        velinds = vel < Fp.maxvelocity;
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