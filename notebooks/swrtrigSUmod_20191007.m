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
%}
pconf = paramconfig;
respwin = [0 250]; % response period in ms rel to swr on
basewin = [-350 -100]; % baseline period in ms rel to swr on
minLBSwr = 10;
minLBSwrSpikes = 10;
numshuffs = 1000;

create_filter = 0;
run_ff = 0;
save_ffdata = 0;
load_ffdata = 0;
combineEpochs = 0;
saveCombinedEpochs = 0;
loadCombinedEpochs = 0;
calcmod = 0;
plotfigs = 1;

pausefigs = 1;
savefigs = 0;
figname = 'lickBoutSUswrmod';
% conditions = {'lickbouts', 'nolickbouts'};
% for c = 1:length(conditions)
%     clear Fp F
%     condition = conditions{c};
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
%     Fp.filtfunction = 'dfa_lickBoutSpikeCorr';
%     Fp.filtfunction = 'dfa_lickswrcorr';
Fp.filtfunction = 'dfa_riptrigspiking';
Fp.params = {'savefigs', 'wtrackdays', 'valid_ntrodes', Fp.filtfunction, 'nonMU_cells'};
Fp = load_filter_params(Fp);

%% FF
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', Fp.tetfilter, ...
        'excludetime', Fp.timefilter,'iterator',Fp.iterator, 'cells', Fp.cellfilter);
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
% ---------------- combine epochs ---------------------------------------------
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

%% selecting ca1 clusters, lickbout swrs, cacluate swr mod
if calcmod
    modF = struct;
    for ian = 1:length(F)
        modF(ian).data = struct;
        cn = 0;
        %% create new data filter for ca1 clusters, swrs during lickbout intervals
        l.animal = F(ian).animal;
        l.params = {'wtrackdays', 'lickbouts', 'ca1SU', 'swrlickmod'};
        l = load_filter_params(l);
        lickca1F = createfilter('animal', l.animal, 'epochs', l.epochfilter,...
            'excludetime', l.timefilter, 'cells', l.cellfilter);
        LBdtc = [];
        for ie = 1:size(lickca1F.epochs{1},1)
            ieTC = lickca1F.data{1}{ie};
            d = repmat(lickca1F.epochs{1}(ie,1), size(ieTC,1),1);
            LBdtc = [LBdtc; d ieTC];
        end
        LBdata = F(ian).data(ismember(cell2mat({F(ian).data.dtc}'), LBdtc, 'rows'));
        for ic = 1:length(LBdata)
            %% make design mat (lickbout-swrs)
            dtc = LBdata(ic).dtc;
            swrStart = LBdata(ic).eventtags(:,2); % col 2 is swr starttime
            epIdx = find(lickca1F(ian).epochs{1}(:,1) == dtc(1));
            DayexIntervals = [];
            for iep = 1:numel(epIdx)
                DayexIntervals = [DayexIntervals; lickca1F.excludetime{1}{epIdx(iep)}];
            end
            lickBoutSwrMask = ~isExcluded(swrStart, DayexIntervals);
            if sum(lickBoutSwrMask) < minLBSwr % filter num lickbout-swrs
                fprintf('%d %d %d only %d swrs in lickbouts. skipping\n',dtc, ...
                    sum(lickBoutSwrMask));
                continue
            end
            %% meanmod score for this cluster, per condition
            ilickFRpsth = LBdata(ic).instantFR(lickBoutSwrMask,:);
            time = LBdata(ic).time;
            % response
            respIdx = [knnsearch(time', respwin(1)/1000) knnsearch(time', respwin(2)/1000)];
            respmeanperSWR = nanmean(ilickFRpsth(:,respIdx(1):respIdx(2)),2)+1e-6;
            % baseline
            baseIdx = [knnsearch(time', basewin(1)/1000) knnsearch(time', basewin(2)/1000)];
            basemeanperSWR = nanmean(ilickFRpsth(:,baseIdx(1):baseIdx(2)),2)+1e-6;
            % mean pct change from baseline
            pctSwrResp = 100*(respmeanperSWR-basemeanperSWR)./basemeanperSWR;
            mod = mean(pctSwrResp);
            
            %% Shuffle
            numSamples = size(ilickFRpsth,2);
            r = randi(numSamples,size(ilickFRpsth,1),numshuffs); % shuff shift offsets
            shuffmod = nan(numshuffs,1);
            for ish = 1:numshuffs
                ishLickFRpsth = circshift(ilickFRpsth, r(:,ish));
                % response
                respShmeanperSWR = nanmean(ishLickFRpsth(:,respIdx(1):respIdx(2)),2)+1e-6;
                % baseline
                baseShmeanperSWR = nanmean(ishLickFRpsth(:,baseIdx(1):baseIdx(2)),2)+1e-6;
                % mean pct change from baseline
                pctSwrResp = 100*(respShmeanperSWR-baseShmeanperSWR)./baseShmeanperSWR;
                shuffmod(ish,1) = nanmean(pctSwrResp);
            end
            %% Test real-mod shuff-rank
            respAboveShuffPct = 100*(1-(sum(shuffmod > mod)./numshuffs));
            fprintf('%s %d %d %d mod is > %.02fperc of shuffs\n', l.animal, dtc, respAboveShuffPct)
            %% ouput
            cn = cn +1;
            modF(ian).animal = l.animal;
            modF(ian).data(cn).dtc = dtc;
            modF(ian).data(cn).time = time;
            modF(ian).data(cn).psth = LBdata(ic).psth;
            modF(ian).data(cn).instantFR = LBdata(ic).instantFR;
            modF(ian).data(cn).lickBoutSwrMask = lickBoutSwrMask;
            modF(ian).data(cn).epoch_noevents = LBdata(ic).epoch_noevents;
            modF(ian).data(cn).mod = mod;
            modF(ian).data(cn).shuffmod = shuffmod;
            modF(ian).data(cn).respAboveShuffPct = respAboveShuffPct;
        end
    end
end
% save calcmod?

%% % SU raster, psth, mod score, population percentage mod scores
if plotfigs
    for ani = 1:length(modF)
        animal = modF(ani).animal;
        andef = animaldef(animal);
        tetinfo = loaddatastruct(andef{2}, animal, 'tetinfo');
        cellinfo = loaddatastruct(andef{2}, animal, 'cellinfo');
        numclusters = length(modF(ani).data);
        for ic = 1:length(numclusters)
            day = modF(an).data(ic).dtc(1);
            nt = modF(an).data(ic).dtc(2);
            clust = modF(an).data(ic).dtc(3);
            Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
                'savefigs', savefigs);
            time = modF(an).data(1).time;
            s = knnsearch(time', -Pp.win(1));
            e = knnsearch(time', Pp.win(2));
           %% SF1 all swr SU raster left
            sf1 = subaxis(4,2,[1 3 5],Pp.posparams{:});
            sf1.Tag = 'SUraster';
            [xx, yy] = find(modF(an).data(ic).psth(:, s:e)');
            scatter(xx/1000-Pp.win(1),yy, 1, '.k')
            ylabel('ripnum')
            xlim([-Pp.win(1) Pp.win(2)])
            ylim([0 size(modF(an).data(ic).psth,1)])
            epoch_noevents = modF(an).data(1).epoch_noevents;
%             epoch_noevents = epoch_noevents(epoch_noevents>0);
            if length(epoch_noevents) > 1
                line(xlim, repmat(epoch_noevents,2,1))
            end
            xticks([])
            
           %% SF2 mean firing rate
            sf2 = subaxis(4,2,7,Pp.posparams{:});
            sf2.Tag = 'meanFR';
            %             [h, e] = histcounts(xx/1000-1.001, -1:.02:1,'Normalization', 'probability');
            %             area(e(1:end-1)+diff(e)/2,h)
            plot(time(s:e), modF(an).data(1).instantFRmean(s:e), 'k')
            axis tight
            ylabel('meanFR Hz')
            xlabel('time s')
            
            %            %% SF3 SU raster right non-bout
            %             sf1 = subaxis(4,2,[1 3 5],Pp.posparams{:});
            %             sf1.Tag = 'SUraster';
            %             [xx, yy] = find(ppF(an).data(iclust).psth(:, s:e)');
            %             scatter(xx/1000-Pp.win(1),yy, 1, '.k')
            %             ylabel('ripnum')
            %             xlim([-Pp.win(1) Pp.win(2)])
            %             ylim([0 size(ppF(an).data(iclust).psth,1)])
            %             epoch_noevents = ppF(an).data(1).epoch_noevents;
            %             epoch_noevents = epoch_noevents(epoch_noevents>0);
            %             if length(epoch_noevents) > 1
            %                 line(xlim, repmat(epoch_noevents,2,1))
            %             end
            %             xticks([])
            %
            %             %% SF4 mean firing rate right non-bout
            %             sf2 = subaxis(4,2,7,Pp.posparams{:});
            %             sf2.Tag = 'meanFR';
            % %             [h, e] = histcounts(xx/1000-1.001, -1:.02:1,'Normalization', 'probability');
            % %             area(e(1:end-1)+diff(e)/2,h)
            %             plot(time(s:e), ppF(an).data(1).instantFRmean(s:e), 'k')
            %             axis tight
            %             ylabel('meanFR Hz')
            %             xlabel('time s')
            
            if 0
                figure
                % check plot compare firing rate psth to spike psth
                subplot(2,2,1)
                imagesc(flipud(zscore(smoothdata(LBdata(ic).instantFR,2,'loess',500), [], 2)))
                ylabel('ripnum')
                colormap(parula)
                title('zfrate')
                
                subplot(2,2,2)
                %             respNumSpikes = sum(icpsth(:,respIdx(1):respIdx(2)));
                [xx, yy] = find(LBdata(ic).psth');
                scatter(xx,yy,1,'k.')
                ylim([0 size(LBdata(ic).psth,1)])
                xlim([0 2002])
                hold on;
                lx = xx(ismember(yy,find(lickBoutSwrMask)));
                ly = yy(ismember(yy,find(lickBoutSwrMask)));
                %             [xx, yy] = find(LBdata(ic).psth(lickBoutSwrMask,:)'); % the psth of lick bout swrs
                scatter(lx,ly,10,'ro')
                yticks([])
                title('red:lickswrs')
                hold off;
                
                subplot(2,2,3)
                imagesc(flipud(zscore(smoothdata(ilickFRpsth,2,'loess',500), [], 2)))
                title('real')
                
                subplot(2,2,4)
                histogram(shuffmod,100)
                axis tight
                yl = ylim;
                hold on;
                line([meanPctSwrResp meanPctSwrResp], [0, yl(2)], 'color', 'r');
                hold off;
            end
            
            
            %% Ripple line
            line(sf1, [0 0], ylim(sf1), 'linestyle', '--', 'color', 'c')
            line(sf2, [0 0], ylim(sf2), 'linestyle', '--', 'color', 'c')
            %% Super Axis
            epoch4subarea = 2;
            area = tetinfo{day}{epoch}{nt}.area;
            subarea = tetinfo{day}{epoch}{nt}.subarea;
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
end
%%

% % separate lick bout swrs from non lick bout swrs
% % where is an example of code snippet to start and end figure plotting?
% an = 1;
% epochs = [6 2; 6 4];
% % get lick bout intervals
% animal = Fp.animals{an};
% andef = animaldef(animal);
% licksvec = getLickBout(andef{2}, animal, epochs);
% lickintervs =
% swrtimes =
%  = isExcluded(lickintervs, swrtimes); % lick bout swrs
% % do i need a position filter to only include non lick bout swrs that occur
% % at well?
%
% [xx, yy] = find(ppF(an).data(1).psth');
% %%
% subplot(4,1,1:3)
% scatter(xx/1000-1.001,yy, '.')
% ylabel('ripnum')
% line(xlim, repmat(ppF(an).data(1).epoch_noevents,2,1))
% xticks([])
% axis tight
% %%
% subplot(4,1,4)
% [h, e] = histcounts(xx/1000-1.001, -1:.02:1,'Normalization', 'probability');
% % zh = zscore(h);
% area(e(1:end-1)+diff(e)/2,h)
% axis tight
% ylabel('probability')
% xlabel('time s')

