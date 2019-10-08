
%{
- Given an input: animal, day, epoch, timestart, timeend
- plot the datachunk of position, speed, lfp example, swrs, rippwrtrace, (su, mu)
%}
pconf = paramconfig;

% d = {'D10', 6, 2, [27 28], 11, 3, 2410, 2419;...
%     'D10', 6, 2, [27 28], 11, 3, 2798, 2804};
d = {'JZ4', 2, 2, [29 30], 3, 13, 2643, 2650;...
    'JZ4', 2, 2, [29 30], 3, 13, 2664, 2673};
% d = {'JZ1', 4, 4, [23 24], 10, 9, 5394, 5403;...
%     'JZ1', 4, 4, [23 24], 10, 9, 5412, 5419};
loaddata = 1;
processData = 1;
plotfigs = 1;
pausefigs = 0;
savefigs = 1;
savefigas = {'png'};

for e = 1:size(d,1)
    
    animal = d{e,1};
    day = d{e,2};
    epoch = d{e,3};
    ca1nt = d{e,4};
    smecnt = d{e,5};
    dmecNt = d{e,6};
    tstart = d{e,7};
    tend = d{e,8};
    %%
    if loaddata
        ntrodes = [ca1nt, smecnt, dmecNt];
        animdef = animaldef(animal);
        swr = loaddatastruct(animdef{2}, animal, 'ca1rippleskons', day);
        eeg = loadeegstruct(animdef{2}, animal, 'eeg', day, epoch, [1:30]);
        %     ripple = loadeegstruct(animdef{2}, animal, 'ripple', day, epoch);
        pos = loaddatastruct(animdef{2}, animal, 'pos', day);
        linpos = loaddatastruct(animdef{2}, animal, 'linpos', day);
        spikes = loaddatastruct(animdef{2}, animal, 'spikes', day);
        cellinfo = loaddatastruct(animdef{2}, animal, 'cellinfo', day);
        ntinfo = loaddatastruct(animdef{2}, animal, 'tetinfo', day);
        DIO = loaddatastruct(animdef{2}, animal, 'DIO', day);
        task = loaddatastruct(animdef{2}, animal, 'task', day);
        lick = loaddatastruct(animdef{2}, animal, 'lick', day);
    end
    %%
    if processData
        % get colors
        arTag = cellfun(@(x) ntinfo{day}{epoch}{x}.area,num2cell([1:30]),'un',0)';
        sbTag = cellfun(@(x) ntinfo{day}{epoch}{x}.subarea,num2cell([1:30]),'un',0)';
        % other params
        maxvelocity = 4;
        validnt = find(cell2mat(cellfun(@(x) ~isempty(x), eeg{day}{epoch}, 'un', 0)));
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
        % Process MU into mustack (ntrode, cluster, spiketime)
        c1Nt = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''ca1'')');
        smNt = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''mec'') && (isequal($subarea, ''supf'') || isequal($subarea, ''nsupf''))');
        dmNt = evaluatefilter(ntinfo{day}{epoch}, 'isequal($valid,''yes'') && isequal($area, ''mec'') && (isequal($subarea, ''deep'') || isequal($subarea, ''ndeep''))');
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
        
        % ripkonspwr
        zrippwr = ((swr{day}{epoch}{1}.powertrace-nanmean(swr{day}{epoch}{1}.powertrace))...
            ./nanstd(swr{day}{epoch}{1}.powertrace))';
        % lfp time
        epstartTime = eeg{day}{epoch}{ntrodes(1)}.starttime;
        srate = eeg{day}{epoch}{ntrodes(1)}.samprate;
        nsamps = length(eeg{day}{epoch}{ntrodes(1)}.data);
        lfptime = [(0:nsamps-1)/srate+epstartTime]'; %repmat([,1,length(ntrodes));
        
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
        
    end
    if plotfigs
        %%
        Pp=load_plotting_params({'defaults','dfa_plotDataChunks', 'examples'}, 'pausefigs', pausefigs, ...
            'savefigs', savefigs);
        %% plot velocity, pos, hd
        subaxis(Pp.nrows,1,1,Pp.posparams{:});
        sfa.Tag = 'top';
        velcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'vel-loess'));
        xcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'x-loess'));
        ycol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'y-loess'));
        dircol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'dir-loess'));
        posIdx = knnsearch(postime(:,1),[tstart; tend]);
        vel = pos{day}{epoch}.data(posIdx(1):posIdx(2),velcol);
        veltime = postime(posIdx(1):posIdx(2));
        %     plot(veltime, vel);
        plot(veltime,vel,'LineWidth',2,'Color', 'k', 'linestyle', '-');
        hold on;
        axis tight; % tighten to get appropriate ylims for stems
        yl = ylim;
        ylabel('vel cm/s');
        set(gca,'YGrid','off','XGrid','off','TickDir','out','TickLength', [0.001 0], ...
            'Xticklabel',[]);
        xlim([tstart tend]);
        xl = xlim;
        xticks(ceil(xl(1)):2:floor(xl(2)))
        %     velinds = vel < maxvelocity;
        %     v = vel(velinds);
        %     v(v<maxvelocity) = 1;
        %     v(v>=maxvelocity) = 0;
        %     stem(veltime(velinds), v*yl(2),'k','filled','Marker', 'none', 'linewidth',2, ...
        %         'color',[.8 .8 .8 .2]);
        hold off
        %% lin pos
        subaxis(Pp.nrows,1,2);
        lpIdx = knnsearch(linpostime(:,1),[tstart; tend]);
        lpClr = linpos{day}{epoch}.statematrix.segmentIndex(lpIdx(1):lpIdx(2));
        linpostimeCh = linpostime(lpIdx(1):lpIdx(2));
        lindistCh = lindist1D(lpIdx(1):lpIdx(2));
        scatter(linpostimeCh,lindistCh,10,'k','filled', 'marker', 'd');
        hold on;
        axis tight;
        %     colormap(Pp.linposcmap);
        ylabel('lin pos cm');
        xlim([tstart tend]);
        xl = xlim;
        %     xticks(ceil(xl(1)):2:floor(xl(2)))
        set(gca,'TickDir','out','YGrid','off', 'XGrid', 'off','TickLength', [0.001 0], ...
            'Tag', 'linpos', 'Xticklabel',[]); %
        %     for iseg = 1:length(segSE(:,1))
        %         patch(sf2, [xl(1) xl(2) xl(2) xl(1)], ...
        %             [segSE(iseg,1) segSE(iseg,1) segSE(iseg,2) segSE(iseg,2)], iseg,...
        %             'FaceAlpha',.05,'edgecolor','none');
        %     end
        %     scatter(linpostimeCh,lindistCh,10,lpClr,'filled', 'marker', 'd', 'markeredgecolor','k');
        hold off;
        
        %% Plot MU spikes
        if any(mustack)
            muCh = mustack(all([mustack(:,3)>tstart mustack(:,3)<tend],2),:);
            if ~isempty(muCh)
                subaxis(Pp.nrows,1,3);
                set(gca, 'Tag', 'mu');
                uclr = zeros(size(muCh,1),3);
                uclr(ismember(muCh(:,1),c1Nt),:)=repmat(Pp.c1clr,sum(ismember(muCh(:,1),c1Nt)),1);
                uclr(ismember(muCh(:,1),dmNt),:)=repmat(Pp.dmclr,sum(ismember(muCh(:,1),dmNt)),1);
                uclr(ismember(muCh(:,1),smNt),:)=repmat(Pp.smclr,sum(ismember(muCh(:,1),smNt)),1);
                dx = scatter(muCh(:,3), muCh(:,4),2,uclr,'filled');
                ylabel('MU');
                set(gca,'TickDir','out','YGrid','off','XGrid','on','fontsize', 10,...
                    'TickLength', [0.001 0],'Tag', 'ndata', 'Xticklabel',[], 'Yticklabel',[]);
                axis tight;
                yl = ylim;
                ylim([yl(1)-2 yl(2)+2]); %
                xlim([tstart tend]);
                xl = xlim;
                %             xticks(ceil(xl(1)):1:floor(xl(2)))
                %             [yt, ydx] = unique(muCh(:,4));
                %             yticks(yt);
                %             yticklabels(muCh(ydx,1));
            end; end
        
        
        %% plot wide lfp of selected ntrodes
        subaxis(Pp.nrows,1,4);
        set(gca, 'Tag', 'lfp');
        lfpIdx = knnsearch(lfptime,[tstart; tend]);
        lfpCh = cell2mat(cellfun(@(x) x.data(lfpIdx(1):lfpIdx(2))', ...
            eeg{day}{epoch}(ntrodes),'un',0)');
        gridmat = Pp.Yoffset*repmat((1:length(lfpCh(:,1)))', 1,length(lfpCh(1,:)));
        eeggrid = (lfpCh + gridmat)';
        plfp = plot(lfptime(lfpIdx(1):lfpIdx(2))', eeggrid,'LineWidth',1);
        ylabel('LFP');
        hold on;
        axis tight;
        yl = ylim;
        xlim([tstart tend]);
        xl = xlim;
        set(gca,'TickDir','out');
        icolors = colorPicker({'mec', 'mec', 'ca1', 'ca1'},...
            'subtags', {'deep', 'supf', 'd', 'd'});
        set(plfp, {'color'}, flipud(icolors'));
        set(gca, 'YGrid', 'off', 'XGrid', 'off','TickLength', [0.001 0], 'Xticklabel',[]);
        ntrodeOffsets = gridmat(:,1);
        yticklabels({'CA1', 'CA1', 'MECs', 'MECd'}); %fliplr(ntrodes));
        yticks([ntrodeOffsets]);
        set(gca,'fontsize',10);
        hold off;
        
        %% plot ripple kons trace
        subaxis(Pp.nrows,1,5);
        set(gca, 'Tag', 'ripkonspwr');
        time = swr{day}{epoch}{1}.eegtimesvec_ref';
        pwrIdx = knnsearch(time,[tstart; tend]);
        % lfp ripple zcored pwr
        rippwrinWin = zrippwr(pwrIdx(1):pwrIdx(2));
        %         noiseinWin = zNoiseRipRatio(pwrIdx(1):pwrIdx(2));
        %         plot(time(pwrIdx(1):pwrIdx(2))',noiseinWin,'LineWidth',2,'Color',[.8 .3 .3 .5]);
        plot(time(pwrIdx(1):pwrIdx(2))', rippwrinWin,'LineWidth',2,'Color',[.3 .3 .3]);
        ylabel('rip pwr z')
        xticks(ceil(xl(1)):2:floor(xl(2)))
        xlabel('Time s');
        axis tight
        %% plot ripple patches
        eventTime = [swr{day}{epoch}{1}.starttime swr{day}{epoch}{1}.endtime];
        swrInWinIdx = all([eventTime(:,1)>tstart eventTime(:,2)<tend],2);
        %     if ~any(swrInWinIdx); fprintf('no rips in win %1.f : %1.f\n', tstart, tend);
        %         if skipNoRipChunks; fprintf('skipping\n');
        %             continue;
        %         end
        %     end
        if any(swrInWinIdx)
            swrInWin = eventTime(swrInWinIdx,:);
            xs = swrInWin(:,1);
            xe = swrInWin(:,2);
            ax = get(gcf, 'children');
            for ix = 1:numel(ax)
                %             if strcmp(ax(ix).Tag, 'image')
                %                 continue;
                %             end
                try
                    %                 disp(ax(ix).Tag)
                    yl = ylim(ax(ix));
                    %                 if strcmp(ax(ix).Tag, 'zrippwr')
                    %                     hold(ax(ix),'on')
                    %                     s = scatter(ax(ix), xs,repmat(yl(2),size(xs,1),1),30,'filled',...
                    %                         'marker', 'd', 'markerfacecolor', 'y', 'markeredgecolor', 'k');
                    %                     hold(ax(ix),'off')
                    %                 end
                    patch(ax(ix), [xs'; xe'; xe'; xs'],repmat([yl(1) yl(1) yl(2) yl(2)]',1,...
                        length(xe)),Pp.ripclr, 'FaceAlpha',Pp.ripFAlpha, 'edgecolor','none');
                catch
                    continue
                end; end; end
        
        %% plot dio stems (licks, rewards)
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
            ax = get(gcf, 'children');
            for ix = 1:length(ax)
                if strcmp(ax(ix).Tag, 'image')
                    continue; end
                try
                    %                 disp(ax(ix).Tag)
                    if strcmp(ax(ix).Tag, 'top')
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
                    l = line(ax(ix),dDF(:,[2 2])',ylim(ax(ix)), 'linewidth',1);
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
        %% ----- link x axis -----
        %     allAxesInFigure = findall(gcf,'type','axes');
        %     linkaxes(allAxesInFigure, 'x');
        % ---- super axis -----
        sprtit = sprintf('lickswrExample %s %d %d %1.f-%1.f', animal,day,epoch,tstart,tend);
        sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', gcf);
        iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
        set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment','center','FontSize',12);
        h = get(gcf,'Children');
        set(gcf,'Children',flip(h)); % put super axis at bottom of axis stack. allows for zoom
        % ---- pause, save figs ----
        if pausefigs; pause; end
        outdir = sprintf('%s/lickswrExample/',pconf.andef{4});
        if savefigs; save_figure(outdir, sprtit, 'savefigas', savefigas); end
    end
end