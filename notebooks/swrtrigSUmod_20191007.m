%{
 SWR-trig SU modulation for 2 conditions:
   - SWRs within lick bout intervals
   - SWRs at well but more than X s away from lick vout intervals

singlecell anal
- 2 timefilters
- how did i do the lick interval thing already.. i think i've saved lick
events so now i can pass them in as a datatype.. in order to get the lick
bout intervals i think i call getLickBout.m for a given epoch..
- i prob don't want to have to do that for all single cells since all
within each epoch will use the same lick bout intervals.. so either i could
save the intervals as another datatype and use singlecellanal..
- or use single epoch anal and nest the single cell loop myself within the
dffunction
- /home/droumis/Src/Matlab/filterframework_dr/notebooks/boutSpikeLFPLocking_20190921.m
- that was how i attempted to do the two timefilter conditions

- fuck... since i want to combine across epochs.. i need to just use the ff
to collect the swrtrig SU spiking, and then in the script i need to compute
mod, plot figs, etc

Text from Neuron paper about SU mod scoring:
SWR-modulation metric: for each PFC neuron we first averaged its SWRtriggered
raster to yield a peri-SWR time histogram (SWR-PSTH). We then created 5000
shuffled SWR-PSTHs, each one constructed by circularly jittering the spikes around
each SWR by a random amount (all spikes around the same SWR were shifted by the
same fixed time). We calculated the squared difference between the real SWR-PSTH
and the mean of the shuffled SWR-PSTHs in the 0-200 ms window after SWR onset, to
obtain a SWR-modulation measure. To determine SWR-modulation significance, we
calculated SWR-modulation measures in the same way for the 5000 shuffled SWRPSTHs.
If the SWR modulation measure of the real SWR-PSTH was greater than 95%
of the shuffled PSTHs, the PFC neuron was determined to be SWR-modulated.
... SWR-modulated neurons were categorized as SWR-excited
or SWR-inhibited by comparing the rate in the 0-200 ms window after SWR onset with
the rate in a pre-SWR background window -500 to -100 ms window before SWR onset.

Text from Mari's paper about SU mod scoring:
To detect significant SWR-modulation of NAc cells, we followed a procedure described
previously.  Briefly, for each cell, we circularly shuffled each SWR-triggered spike
train by a random amount up to ±0.5 s to generate 5000 shuffled PETHs.  We then
calculated the summed squared difference of the real PETH relative to the mean of
the shuffles in a 0-200 ms window post SWR-onset, and compared it to the same value
 for each shuffle relative to the mean of the shuffles.  Significance at p<0.05
indicates that the real modulation exceeded 95% of the shuffles.  The direction of
 modulation was defined from a modulation index, calculated as the mean firing rate
in the 0-200 ms window minus the mean baseline firing rate from -500 to -100 ms, divided
by the mean baseline firing rate.  This sign of this index was used to assign cells as
significantly positively or negatively SWR-modulated.

kenny doesn't have any real description of his modulation testing.. fucking
idiot/

gideon used a different response period.. -250:250. and did 1k shuffles
 
for each cluster:
1. circ shuffle rand 0:500ms swr psth, save 0-200 mean rate for each shuffle
2. squared diff of real 0-200ms mean vs squared shuffle 0-200 mean of means
3. is #2 value > 95% of 0-200ms shuffle means? ~ significantly modulated
4. sign of real 0-200ms mean rate minus real -500:-100ms determines excited
or inhibited, if also significantly modulated via #3..

!! because this is continuous rate.. i can treat it like lfp.. and do that
perm test in a similar way to test the spectral diff of event trig response
between two conditions.. i would only need to use like 2:50 Hz range..

i like this because it will test for differences between conditions, and
against noise, in time-frequency space.. instead of just doing a single
time point as the 'modulation window'..
- what wavelet script was i using to do all this perm testing?
- should i just do the normal way first? that's probably the right thing to
do to verify the % of modulated ca1 neurons vs the other papers.
- what was the % results for this in the other papers?

.. actually this is very similar to the wavelet lfp analysis stuff i was
doing because there i also had groups of indices corresponding to each
condition
First, calculate the full mod score for a given cluster using all the swrs against noise
Next, separate the swr psth's into lick and nonlick swr psth's
    - do this by testing the starttime of each swr with inclusion in lick
    bout intervals
Finally, calculate the modulation score for each condition

finally, finally, make a heatraster of meanFR for all clusters per area (ca1, mecd, mecs)
    - for each condition (lick, nonlick)
.. putting the lick, nonlick, allswrs together, side by side, is probably
    most of the way towards arguing that these are real swr's.. or if they
    are different, that's interesting, but i'll prob need to do decoding

work continued in /home/droumis/Src/Matlab/filterframework_dr/notebooks/sumod_20191013.m
%}
%%
pconf = paramconfig;
respwin = [0 200]; % response period in ms rel to swr on
basewin = [-300 -100]; % baseline period in ms rel to swr on
minNumSwr = 10;
% minNumSwrSpikes = 10;
nshuffs = 1000;
shuffms = 700;
dmatIdx = {'all', 'lickbout'};

create_filter = 0;
run_ff = 0;
save_ffdata = 0;
load_ffdata = 0;
combineEpochs = 0;
saveCombinedEpochs = 0;
loadCombinedEpochs = 0;
calcmod = 0;
loadmodF = 0;
plotClusterFigs = 0;
plotPopulationFigs = 1;

pausefigs = 1;
savefigs = 0;
figname = 'lickBoutSUswrmod';

Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_riptrigspiking';
Fp.params = {'savefigs', 'wtrackdays', 'valid_ntrodes', Fp.filtfunction, 'nonMU_cells'};
Fp = load_filter_params(Fp);

%% FF
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes',...
        Fp.tetfilter, 'excludetime', Fp.timefilter,'iterator',Fp.iterator, 'cells',...
        Fp.cellfilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff; F = runfilter(F);
    for d = 1:length(F); F(d).datafilter_params = Fp;
    end
end
if save_ffdata
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', sprintf('_%s_%s', Fp.epochEnvironment))
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', sprintf('_%s', Fp.epochEnvironment));
end
%% ---------------- combine cluster epochs --------------------------------------
if combineEpochs
    if exist('F', 'var')
        ppF = combine_epochs(F, Fp, saveCombinedEpochs, Fp.paths);
    else
        error('create or load data filter output to combine epochs \n')
    end
end
% ---------------- Load combined epochs ---------------------------------------
if loadCombinedEpochs
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals,...
        'filetail', '_combEps');
end
%% ---------------- calc mod---------------------------------------
if calcmod
    % create new modF data filtered for certain clusters, swr-designmat
    filtF.animal = {F.animal};
    filtF.params = {'wtrackdays', 'lickbouts', 'ca1SU', 'swrlickmod'};
    filtF = load_filter_params(filtF);
    filtF = createfilter('animal', filtF.animal, 'epochs', filtF.epochfilter,...
        'excludetime', filtF.timefilter, 'cells', filtF.cellfilter);
    
    modF = calcSUmod(F,filtF, 'respwin', respwin, 'basewin', basewin, 'minNumSwr', ...
        minNumSwr, 'nshuffs', nshuffs, 'shuffms', shuffms, ...
        'filetail', ['_' Fp.epochEnvironment]);
end
% ---------------- load mod---------------------------------------
if loadmodF
    modF = load_data([pconf.andef{3} '/sumod'], 'sumod', Fp.animals,...
        'filetail', ['_' Fp.epochEnvironment]);
end

%% % SU raster, psth, mod score, population percentage mod scores
% for each cluster,
% dmat conditions:
% 1. plot : all swr SPraster, SPpsth, FRheatrast, FRpsth, m%chVshuffdistr
% 2. plot : lB  swr SPraster, SPpsth, FRheatrast, FRpsth, m%chVshuffdistr
% for each, annotate:
%  : adtc, env, backidx, respidx, riptime, m%change, m%chShuffRank

% All animal ca1 population figs:
% 1. heatraster ranked by m%chShuffRank of smooth zcore psthFR cluster x time
% 2. cdf m%change

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
if plotPopulationFigs
    % for each animal, gather the ca1 su meanFRpsth, sort by score, and
    % heatmap
    for ian = 1:length(modF)
        animal = modF(ian).animal;
        Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
            'savefigs', savefigs);
        time = modF(ian).data(1).time;
        s = knnsearch(time', -Pp.win(1));
        e = knnsearch(time', Pp.win(2));
        %%
        id = 1;
        for ic = 1:length(modF(ian).data)
            FRrast(ic,:) = nanmean(modF(ian).data(ic).inputdata.instantFR(...
                modF(ian).data(ic).dmat(:,id),s:e));
            clustRank(ic,1) = modF(ian).data(ic).respAboveShuffPct{id};
        end
        frsm = smoothdata(zscore(FRrast,[],2),2,'loess',100);
        [~, sortidx] = sort(clustRank, 'descend');
        imagesc(frsm(sortidx,:))
        colormap(jet)
        %% now plot cum dist scores with sig lines a la sirota
        
        
        %% Super Axis
        epoch4subarea = 2;
        area = tetinfo{day}{epoch4subarea}{nt}.area;
        subarea = tetinfo{day}{epoch4subarea}{nt}.subarea;
        %             meanrate = cellinfo{day}{epoch}{nt}{clust}.firingrate;
        %             isolation = cellinfo{day}{epoch}{nt}{clust}.isolation;
        sprtit = setSuperAxTitle(animal, figname, 'addtitle', ...
            sprintf('%d %d %d %s %s', day, nt, clust, area, subarea));
        if pausefigs
            pause;
        end
        if savefigs
            outdir = sprintf('%s/%s/', pconf.andef{4}, figname);
            save_figure(outdir, sprtit);
        end
    end
    
end

