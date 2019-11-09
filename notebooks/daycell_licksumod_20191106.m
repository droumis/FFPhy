
% get data
pconf = paramconfig;
create_filter = 0;
run_ff = 0;
load_ffdata = 0;

% plot
plotPerCell = 1;
pausefigs = 1;
savefigs = 0;


%% filter params
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_lickXCorrSpikes';
Fp.Label = 'lickSpikeXcorr';
Fp.params = {'wtrack', '<4cm/s', 'wtrackdays', 'valid_ntrodes', 'nonMU_cells', Fp.filtfunction};
Fp = load_filter_params(Fp);

%% FF
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes',...
        Fp.tetfilter, 'excludetime', Fp.timefilter,'iterator',Fp.iterator, 'cells',...
        Fp.cellfilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F), 'un', 1);
    F = runfilter(F);
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.env '_' Fp.eventType]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', ['_' Fp.env '_' Fp.eventType]);
end

%% plot
if plotPerCell
    % i want to later include spatial FRmap, cellinfo, swrmod
        % TimeFRheatrast+ModScore/Sig      PhaseFRheatrast
        % Xcorr vs Shuf                    PSPH vs Shuf
        % Excorr vs Shuf                   LogMagMRV vs Shuf
    for a = 1:length(F) % animal
        animal = F(a).animal{3};
        for ic = 1:length(F(a).output{1}) % cell
            Pp=load_plotting_params({'defaults','dfa_lickXCorrSpikes'});
            ncols = 2; 
            nrows = 3;
            day = F(a).output{1}(ic).index(1);
            eps = F(a).output{1}(ic).index(4:5);
            nt = F(a).output{1}(ic).index(2);
            clust = F(a).output{1}(ic).index(3);
            fprintf('%s %d %d %d\n', animal, day, nt, clust);
            time = F(a).output{1}(ic).time';
            s = knnsearch(time, -Pp.win(1));
            e = knnsearch(time, Pp.win(2));
            time = time(s:e);
            %% Time FRheatrast
            sf = subaxis(nrows,ncols, 1, Pp.posparams{:});
            ifrtime = F(a).output{1}(ic).evTrigSpike.time;
            frs = knnsearch(ifrtime', -Pp.win(1));
            fre = knnsearch(ifrtime', Pp.win(2));
            ifrtime = ifrtime(frs:fre);
            try
                FRrast = flipud(F(a).output{1}(ic).evTrigSpike.instantFR(:,frs:fre));
            catch
                fprintf('evTrigSpike is empty \n');
                continue
            end
            FRrastproc = zscore(smoothdata(smoothdata(FRrast,2,'loess',250),1,'loess',10), [], 2);
            imagesc(ifrtime, 1:size(FRrastproc,1), FRrastproc)
            ylabel('lick num')
            colormap(gray)
            line([0 0], ylim, 'linestyle', '-', 'color', [1 .5 .5 .5])
            %% Phase FRheatrast
            sf = subaxis(nrows,ncols, 2, Pp.posparams{:});
            
            %%  Xcorr
            sf = subaxis(nrows,ncols, 3, Pp.posparams{:});
            sf.Tag = 'xcorr';
            xc = F(a).output{1}(ic).smthxc(s:e)';
            p1 = plot(time, xc, 'b', 'linewidth', 1, 'Displayname', 'xcorr');
            p1.Color = [0 0 0 0];
            hold on
            axis tight
            yl = ylim;
            xcSh = cell2mat(F(a).output{1}(ic).smthxcShuf');
            xcSh = xcSh(:,s:e);
            xcShM = nanmean(xcSh)';
            xcShSTD = nanstd(xcSh)';
            a1 = area(time, xc, (yl(2)-yl(1))/1.5);
            a1.FaceColor = 'b';
            a1.FaceAlpha = .5;
            a1.EdgeAlpha = 0;
            plot(time, xcShM, 'color', 'k', 'linewidth', 1); hold on;
            fill([time; flipud(time)],[xcShM-xcShSTD;flipud(xcShM+xcShSTD)],'k',...
                'linestyle','none', 'facealpha', .3);
            line([0 0],ylim, 'color','red', 'linewidth',2, 'Color', [1 0 0 .5]);
            hold off
            ylabel('xcorr')
            xlabel('lag time s')
            set(gca, 'YGrid', 'on', 'XGrid', 'on','TickDir','out', 'TickLength', [0.001 0]);
            
            %% PeriLickPhaseHistogram vs Shuf
            sf = subaxis(nrows,ncols, 4, Pp.posparams{:});
            
            %% Excorr vs Shuf
            sf = subaxis(nrows,ncols, 5, Pp.posparams{:});
            exc = F(a).output{1}(ic).excorr';
            excSh = cell2mat(F(a).output{1}(ic).excorrShuf');
            histogram(excSh, 60,'Normalization','probability','edgealpha', 0, 'facecolor', 'k');
            hold on;
            excShsort = sort(excSh);
            idxsig = round(Fp.sigpct/100*length(excShsort));
            line([excShsort(idxsig) excShsort(idxsig)], ylim, 'color', [0 0 0 .8], 'linestyle', '--');
            idxsig = round((1-(Fp.sigpct/100))*length(excShsort));
            line([excShsort(idxsig) excShsort(idxsig)], ylim, 'color', [0 0 0 .8], 'linestyle', '--');
            line([exc exc], ylim, 'color', 'r');
            excp = 1-sum(exc>excSh)/length(excSh);
            title(sprintf('excorr %.03f p%.03f', exc, excp));
            ylabel('probability')
            xlabel('excess corr')
            axis tight
            hold off;
            %% LogMagMRV vs Shuf
            
            
            
            %% ---- pause, save figs ----
            if pausefigs; pause; end
            if saveplots
                [~, fname,~] = fileparts(mfilename('fullpath'));
                outdir = sprintf('%s/%s/%s/', pconf.andef{4},fname,animal);
                save_figure(outdir, sprtit, 'savefigas', savefigas);
            end
        end
    end
end


if plotClusterFigs
    for ian = 1:length(modF)
        animal = modF(ian).animal;
        andef = animaldef(animal);
        tetinfo = loaddatastruct(andef{2}, animal, 'tetinfo');
        cellinfo = loaddatastruct(andef{2}, animal, 'cellinfo');
        numclusters = length(modF(ian).data);
        for ic = 1:numclusters
            day = modF(ian).data(ic).dtc(1);
            nt = modF(ian).data(ic).dtc(2);
            clust = modF(ian).data(ic).dtc(3);
            %             if any(cell2mat(modF(ian).data(ic).numSpikesResp)) < 10
            %                 fprintf('%d %d %d cond %d not enough spikes in response period\n', ...
            %                     day,nt,clust, id);
            %                 continue
            %             end
            %             if any(cell2mat(modF(ian).data(ic).numSpikesBase)) < 10
            %                 fprintf('%d %d %d cond %d not enough spikes in baseline period\n', ...
            %                     day,nt,clust);
            %                 continue
            %             end
            Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
                'savefigs', savefigs);
            time = modF(ian).data(1).time;
            s = knnsearch(time', -Pp.win(1));
            e = knnsearch(time', Pp.win(2));
            %             for id = 1:size(modF(ian).data(ic).dmat,2)
            %% swrtrig SU spike raster
            id = 1;
            sf1 = subaxis(6,3,[1 4],Pp.posparams{:});
            sf1.Tag = 'SUraster';
            [xx, yy] = find(modF(ian).data(ic).inputdata.psth(find(modF(ian).data(ic).dmat(:,id)), s:e)');
            scatter(xx/1000-Pp.win(1),yy, 10, '.k', 'markeredgealpha', .6)
            hold on;
            ylabel('swr')
            xlim([-Pp.win(1) Pp.win(2)])
            ylim([0 size(modF(ian).data(ic).inputdata.psth(modF(ian).data(ic).dmat(:,id),:),1)])
            %                 epoch_noevents = modF(ian).data(ic).inputdata.epoch_noevents;
            %                 if length(epoch_noevents) > 1
            %                     line(xlim, repmat(epoch_noevents,2,1))
            %                 end
            xticks([])
            yl = ylim;
            patch([0 .2 .2 0], [yl(1) yl(1) yl(2) yl(2)], [.2 .6 .2], 'facealpha', .1, 'edgealpha', 0)
            patch([-.3 -.1 -.1 -.3], [yl(1) yl(1) yl(2) yl(2)], [.6 .2 .9], 'facealpha', .1, 'edgealpha', 0)
            line([0 0], ylim, 'linestyle', '--', 'color', [1 .5 .5 .5])
            %                 line([pctConfLow pctConfLow], ylim, 'linestyle', ':', 'color', [.5 .5 .5 .5])
            hold off
            
            % spike psth
            sf2 = subaxis(6,3,7,Pp.posparams{:});
            sf2.Tag = 'psth';
            bintime = linspace(time(s), time(e), 50);
            %                 [h, ~] = histcounts(xx/1000-Pp.win(1),bintime,'Normalization',...
            %                     'probability');
            histogram(xx/1000-Pp.win(1),bintime,'Normalization',...
                'probability', 'edgealpha', 0, 'facecolor', 'k', 'facealpha', 1)
            %                 binc = bintime(1:end-1)+diff(bintime(1:2));
            %                 area(binc,h)
            axis tight
            hold on;
            yl = ylim;
            xticks([])
            patch([0 .2 .2 0], [yl(1) yl(1) yl(2) yl(2)], [.2 .6 .2], 'facealpha', .1, 'edgealpha', 0)
            patch([-.3 -.1 -.1 -.3], [yl(1) yl(1) yl(2) yl(2)], [.6 .2 .9], 'facealpha', .1, 'edgealpha', 0)
            line([0 0], ylim, 'linestyle', '--', 'color', [1 .5 .5 .5])
            hold off;
            
            %% lickbout swrtrig SU spike raster
            id = 2;
            sf1 = subaxis(6,3,[10 13],Pp.posparams{:});
            sf1.Tag = 'SUraster';
            [xx, yy] = find(modF(ian).data(ic).inputdata.psth(find(modF(ian).data(ic).dmat(:,id)), s:e)');
            scatter(xx/1000-Pp.win(1),yy, 10, '.k', 'markeredgealpha', .6)
            hold on;
            ylabel('swr')
            xlim([-Pp.win(1) Pp.win(2)])
            ylim([0 size(modF(ian).data(ic).inputdata.psth(modF(ian).data(ic).dmat(:,id),:),1)])
            %                 epoch_noevents = modF(ian).data(ic).inputdata.epoch_noevents;
            %                 if length(epoch_noevents) > 1
            %                     line(xlim, repmat(epoch_noevents,2,1))
            %                 end
            xticks([])
            yl = ylim;
            patch([0 .2 .2 0], [yl(1) yl(1) yl(2) yl(2)], [.2 .6 .2], 'facealpha', .1, 'edgealpha', 0)
            patch([-.3 -.1 -.1 -.3], [yl(1) yl(1) yl(2) yl(2)], [.6 .2 .9], 'facealpha', .1, 'edgealpha', 0)
            %                 line([pctConfLow pctConfLow], ylim, 'linestyle', ':', 'color', [.5 .5 .5 .5])
            line([0 0], ylim, 'linestyle', '--', 'color', [1 .5 .5 .5])
            hold off
            
            % spike psth
            sf2 = subaxis(6,3,16,Pp.posparams{:});
            sf2.Tag = 'psth';
            bintime = linspace(time(s), time(e), 50);
            %                 [h, ~] = histcounts(xx/1000-Pp.win(1),bintime,'Normalization',...
            %                     'probability');
            histogram(xx/1000-Pp.win(1),bintime,'Normalization',...
                'probability', 'edgealpha', 0, 'facecolor', 'k', 'facealpha', 1)
            %                 binc = bintime(1:end-1)+diff(bintime(1:2));
            %                 area(binc,h)
            axis tight
            hold on;
            yl = ylim;
            patch([0 .2 .2 0], [yl(1) yl(1) yl(2) yl(2)], [.2 .6 .2], 'facealpha', .1, 'edgealpha', 0)
            patch([-.3 -.1 -.1 -.3], [yl(1) yl(1) yl(2) yl(2)], [.6 .2 .9], 'facealpha', .1, 'edgealpha', 0)
            line([0 0], ylim, 'linestyle', '--', 'color', [1 .5 .5 .5])
            hold off;
            
            %% FR heatraster
            id = 1;
            sf3 = subaxis(6,3,[2 5],Pp.posparams{:});
            sf3.Tag = 'FRheatraster';
            FRrast = flipud(modF(ian).data(ic).inputdata.instantFR(modF(ian).data(ic).dmat(:,id),s:e));
            FRrastproc = zscore(smoothdata(smoothdata(FRrast,2,'loess',250),1,'loess',10), [], 2);
            imagesc(time(s:e), 1:size(FRrastproc,1), FRrastproc)
            %                 ylabel('ripnum')
            colormap(jet)
            mod = modF(ian).data(ic).meanPctSwrResp{id};
            title(sprintf('%.02f pctbaseline', mod))
            yticks([])
            xticks([])
            line([0 0], ylim, 'linestyle', '--', 'color', [1 .5 .5 .5])
            
            % mean FR
            sf4 = subaxis(6,3,8,Pp.posparams{:});
            sf4.Tag = 'meanFR';
            area(time(s:e), smoothdata(nanmean(FRrast),2,'loess',250), 'facecolor', 'k')
            % add conf bound fills
            axis tight
            %                 ylabel('meanFR Hz')
            yticks([])
            xticks([])
            hold on;
            yl = ylim;
            patch([0 .2 .2 0], [yl(1) yl(1) yl(2) yl(2)], [.2 .6 .2], 'facealpha', .1, 'edgealpha', 0)
            patch([-.3 -.1 -.1 -.3], [yl(1) yl(1) yl(2) yl(2)], [.6 .2 .9], 'facealpha', .1, 'edgealpha', 0)
            line([0 0], ylim, 'linestyle', '--', 'color', [1 .5 .5 .5])
            hold off;
            
            %% lickbout FR heatraster
            id = 2;
            sf3 = subaxis(6,3,[11 14],Pp.posparams{:});
            sf3.Tag = 'FRheatraster';
            FRrast = flipud(modF(ian).data(ic).inputdata.instantFR(...
                modF(ian).data(ic).dmat(:,id),s:e));
            FRrastproc = zscore(smoothdata(smoothdata(FRrast,2,'loess',250),1,...
                'loess',10), [], 2);
            imagesc(time(s:e), 1:size(FRrastproc,1), FRrastproc)
            %                 ylabel('ripnum')
            colormap(jet)
            mod = modF(ian).data(ic).meanPctSwrResp{id};
            title(sprintf('%.02f pctbaseline', mod))
            yticks([])
            xticks([])
            line([0 0], ylim, 'linestyle', '--', 'color', [1 .5 .5 .5])
            
            % mean FR
            sf4 = subaxis(6,3,17,Pp.posparams{:});
            sf4.Tag = 'meanFR';
            area(time(s:e), smoothdata(nanmean(FRrast),2,'loess',250), 'facecolor', 'k')
            % add conf bound fills
            axis tight
            %                 ylabel('meanFR Hz')
            xlabel('time s')
            yticks([])
            hold on;
            yl = ylim;
            patch([0 .2 .2 0], [yl(1) yl(1) yl(2) yl(2)], [.2 .6 .2], 'facealpha', .1, 'edgealpha', 0)
            patch([-.3 -.1 -.1 -.3], [yl(1) yl(1) yl(2) yl(2)], [.6 .2 .9], 'facealpha', .1, 'edgealpha', 0)
            line([0 0], ylim, 'linestyle', '--', 'color', [1 .5 .5 .5])
            hold off;
            
            %% modShuffRnk
            id = 1;
            sf5 = subaxis(6,3,[3 6 9],Pp.posparams{:});
            sf5.Tag = 'modShuffRnk';
            %                 bintime = linspace(time(s), time(e), 50);
            shuffs = modF(ian).data(ic).shuffMeanPctSwrResp{id};
            %                 [h, bins] = histcounts(shuffs, 100, 'Normalization', 'probability');
            histogram(shuffs,50,'Normalization',...
                'probability', 'edgealpha', 0, 'facecolor', 'k')
            %                 area(bins(1:end-1)+diff(bins(1:2)),h)
            %                 axis tight
            %
            % %                 [h, ~] = histcounts(modF(ian).data(ic).shuffMeanPctSwrResp{id}, 10);
            mod = modF(ian).data(ic).meanPctSwrResp{id};
            xl = xlim;
            [sortShuffs, sortIndex] = sort(shuffs);
            pctConfHigh = sortShuffs(floor(0.975 * numel(shuffs)));
            pctConfLow = sortShuffs(floor(0.025 * numel(shuffs)));
            line([pctConfHigh pctConfHigh], ylim, 'linestyle', ':', 'color', [.5 .5 .5 .5])
            line([pctConfLow pctConfLow], ylim, 'linestyle', ':', 'color', [.5 .5 .5 .5])
            if mod < xl(2)
                hold on;
                line([mod mod], ylim, 'linestyle', '--', 'color', [.5 .5 1 .8])
            end
            title(sprintf('mod > %.0fpct of shuff',modF(ian).data(ic).respAboveShuffPct{id}));
            yticks([])
            % why isn't this centered on 0 pcnt?
            %             xlabel('pct change from baseline')
            hold off
            %% lickbout modShuffRnk
            id = 2;
            sf5 = subaxis(6,3,[12 15 18],Pp.posparams{:});
            sf5.Tag = 'modShuffRnk';
            %                 bintime = linspace(time(s), time(e), 50);
            shuffs = modF(ian).data(ic).shuffMeanPctSwrResp{id};
            %                 [h, bins] = histcounts(shuffs, 100, 'Normalization', 'probability');
            histogram(shuffs,50,'Normalization',...
                'probability', 'edgealpha', 0, 'facecolor', 'k')
            %                 area(bins(1:end-1)+diff(bins(1:2)),h)
            %                 axis tight
            %
            % %                 [h, ~] = histcounts(modF(ian).data(ic).shuffMeanPctSwrResp{id}, 10);
            mod = modF(ian).data(ic).meanPctSwrResp{id};
            xl = xlim;
            [sortShuffs, sortIndex] = sort(shuffs);
            pctConfHigh = sortShuffs(floor(0.975 * numel(shuffs)));
            pctConfLow = sortShuffs(floor(0.025 * numel(shuffs)));
            line([pctConfHigh pctConfHigh], ylim, 'linestyle', ':', 'color', [.5 .5 .5 .5])
            line([pctConfLow pctConfLow], ylim, 'linestyle', ':', 'color', [.5 .5 .5 .5])
            if mod < xl(2)
                hold on;
                line([mod mod], ylim, 'linestyle', '--', 'color', [.5 .5 1 .8])
            end
            title(sprintf('mod > %.0fpct of shuff',modF(ian).data(ic).respAboveShuffPct{id}));
            yticks([])
            % why isn't this centered on 0 pcnt?
            xlabel('pct change from baseline')
            hold off
            
            line(sf2, [0 0], ylim(sf2), 'linestyle', '--', 'color', [1 .5 .5 .5])
            line(sf3, [0 0], ylim(sf3), 'linestyle', '--', 'color', [1 .5 .5 .5])
            line(sf4, [0 0], ylim(sf4), 'linestyle', '--', 'color', [1 .5 .5 .5])
            %% Super Axis
            epoch4subarea = 2;
            barea = tetinfo{day}{epoch4subarea}{nt}.area;
            subarea = tetinfo{day}{epoch4subarea}{nt}.subarea;
            %             meanrate = cellinfo{day}{epoch}{nt}{clust}.firingrate;
            %             isolation = cellinfo{day}{epoch}{nt}{clust}.isolation;
            sprtit = setSuperAxTitle(animal, figname, 'addtitle', ...
                sprintf('%d %d %d %s %s %s', day, nt, clust, barea, subarea, ...
                dmatIdx{id}));
            if pausefigs
                pause;
            end
            if savefigs
                outdir = sprintf('%s/%s/', pconf.andef{4}, figname);
                save_figure(outdir, sprtit);
            end
        end
    end
end