%{
NB 20191117
finalize swr-trig and lick trig results all animals combined
.. continuation of suresplickvsswr_20191116
- *** valid_ntrodes filter is not currently being respected for spiking

eventType: lick, ca1rippleskons
eventSet: all, nonburst, lickburst, wet, dry
data: spikes, eeg
area: ca1, mecS, mecD
enb: wtrack

Figure 2. lick & swr triggered neural activity in ca1 and mec

A.
eventTrigSpiking::
    - barn.rat.beer.wheelbarrow
    // dfa_eventTrigSpiking (rat::F.output{1}(perSU)::time, psfr, eventTimes) ->
        -> makeExpvarCatDesignMat (beer::dmat.dayeps, dm) ->
        -> calcSUMod (wheelbarrow::modF.output{1}(perSU)::time, evMean, mPctChange, modPctRank) ->
                - MOD SCORE: per su, eSet: mean % PSFR resp (0:200) from baseline(-300:-100)
                    - SH: per event, shuf time (700ms)
        -> MAKE plotHeatRaster, plotCDFsig
            -  su plotting past: suSWRvsLick_20191112
    
B.
eventTrigLFP::
    - forest.bear.cactus.mushroom.beer.leaf
    // dfa_eventTrigLFP (bear:: ) ->
        -> stack_lfp (cactus:: ) ->
        -> computeAnalyticSignal (mushroom:: ) ->
        -> makeEventSet (beer:: ) ->
        -> get_power (cactus:: eventSetPowerTrace) ->

%}

pconf = paramconfig;
eventTrigLFP = 0; % forest.bear.cactus.mushroom.beer.leaf == eventSet mean Spect
eventTrigSpiking = 1; % barn.rat.beer.wheelbarrow == eventSet SU mod
eventType = 'lick'; %lick swr

% run FF
create_filter = 0;
run_ff = 0;
load_ffdata = 0;

%% create condition design mat
make_expvarCat = 0; % beer
load_expvarCat = 0;

%% spike
calcSUMod = 0; % wheelbarrow
loadSUMod = 1; 
gatherResults = 1; % combines across animals. keeps eventSet results seperate per unit

%% LFP
% stack data
stack_LFP = 0;
load_LFPstack = 0;
% get power
make_rawpwr = 0;
load_rawpwr = 0;
% compute per condition
make_expvarCatMeanPwr = 0;
load_expvarCatMeanPwr = 0;
% combine per area
% combineArea = 0;

%% plot
plotfigs = 1;
showfigs = 1;
pausefigs = 1;
savefigs = 1;

% swr
plotEventTrigSU = 0;         % per eventSet, per SU
% plotEventTrigHeatrastPerAni = 0; % per eventSet, per area, per animal and all animals
plotEventTrigHeatrastAllAni = 1; % requires gatherResults
plotEventTrigModCDF = 1;     % requires gatherResults// per eventSet, per area, per animal and all animals

% lick
plotLickTrigSU = 0;         % per eventSet, per SU
plotLickTrigHeatRaster = 0; % per eventSet, per area, per animal and all animals
plotLickTrigModCDF = 0;     % per eventSet, per area, per animal and all animals

%%
Fp = [];
Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'}; %, };
Fp.areas = {{'ca1', 'd'}, {'mec', 'deep'}, {'mec', 'supf'}};

if eventTrigLFP
    Fp.filtfunction = 'dfa_eventTrigLFP'; % Bellicose Bear
    if strcmp(eventType, 'lick')
        Fp.Label = 'wtrackLickTrigLFP';
    elseif strcmp(eventType, 'swr')
        Fp.Label = 'wtrackSWRTrigLFP';
    end
    Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
        'excludeAfterLastWell', 'referenced', '4-350Hz',  Fp.Label, Fp.filtfunction};
    wp = getWaveParams(Fp.waveSet);
    
elseif eventTrigSpiking
    Fp.filtfunction = 'dfa_eventTrigSpiking'; % Redolent Rat
    if strcmp(eventType, 'lick')
        Fp.Label = 'wtrackLickTrigSpiking';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
            'excludeAfterLastWell', 'nonMU_cells', Fp.Label, Fp.filtfunction};
    elseif strcmp(eventType, 'swr') %'excludeNoise', 
        Fp.Label = 'wtrackSWRTrigSpiking';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
            'excludeAfterLastWell', 'nonMU_cells', 'ripples', ...
            Fp.Label, Fp.filtfunction}; % 'excludeNoise',
    end
end
Fp = load_filter_params(Fp);

%%
if create_filter
    if eventTrigLFP
        F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', ...
            Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    elseif eventTrigSpiking
        F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes',...
            Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator, 'cells',...
            Fp.cellfilter);
    end
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F), 'un', 1);
    F = runfilter(F);
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.Label]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', ['_' Fp.Label]);
end

%% make design mat to slice the rawpwr trials
if make_expvarCat
    outdir = 'expvarCatBurst';
    expvarCat = makeExpvarCatDesignMat(lfpstack, 'outdir', outdir, 'expvars', {'all'}, ...
        'lfptype', Fp.uselfptype, 'eventType', Fp.eventType);
end
if load_expvarCat
    outdir = 'expvarCatBurst';
    outpath = [pconf.andef{2},outdir,'/'];
    expvarCat = load_data(outpath, [outdir,'_',Fp.env,'_',Fp.eventType], Fp.animals);
end

%% calc su modulation
if calcSUMod % wheelbarrow
    modF = calcSUmod(F);
    save_data(modF, 'results', Fp.Label);
end
if loadSUMod
    modF = load_data('results', Fp.Label, Fp.animals);
end

%% Gather all su per area, eventSet
if gatherResults
    FRHeatrast = {};
    allmodF = cell2mat(arrayfun(@(x) modF(x).output{1}', 1:length(modF), 'un', 0)');
    numESet = length(modF(1).dmatIdx);
    for ar = 1:length(Fp.areas) % per area
        % find cells in this area
        areaIdx = strcmp({allmodF.area}', Fp.areas{ar}{1});
        subareaIdx = ~cellfun(@isempty, strfind({allmodF.subarea}', Fp.areas{ar}{2}), 'un', 1);
        iareaIdx = find(all([areaIdx subareaIdx],2));
        iareaIdx = iareaIdx(arrayfun(@(x) ~isempty(allmodF(x).mPctChange), ...
            iareaIdx,'un',1));
        
        for iv = 1:numESet % per eventSet    
            % modPctRank
            iMPRank{ar,iv} = cell2mat(arrayfun(@(x) allmodF(x).modPctRank{iv}, iareaIdx,'un',0));
            iMPctCh{ar,iv} = cell2mat(arrayfun(@(x) allmodF(x).mPctChange{iv}, iareaIdx,'un',0));
            iMPctChSh{ar,iv} = cell2mat(arrayfun(@(x) allmodF(x).mPctChangeSh{iv}, iareaIdx,'un',0));
            iMZCh{ar,iv} = cell2mat(arrayfun(@(x) allmodF(x).mZChange{iv}, iareaIdx,'un',0));
            iMZChSh{ar,iv} = cell2mat(arrayfun(@(x) allmodF(x).mZChangeSh{iv}, iareaIdx,'un',0));
            iMZResp{ar,iv} = cell2mat(arrayfun(@(x) allmodF(x).mZResp{iv}, iareaIdx,'un',0));
            
            % Firing Rate HeatRaster concat, normalize, smooth, sort
            iFRHRz = cell2mat(arrayfun(@(x) allmodF(x).evMeanZ{iv}, iareaIdx,'un',0));
            iFRHRzsm = smoothdata(iFRHRz, 2, 'loess', 10);
            pctChange = cell2mat(arrayfun(@(x) allmodF(x).mPctChange{iv},iareaIdx,'un',0));
            [~, srtIdx] = sort(pctChange, 1, 'descend');
            iFRHRsmzSorted =  iFRHRzsm(srtIdx,:);
            try
                FRHeatrast{ar,iv} = [FRHeatrast{ar,iv}; iFRHRsmzSorted];
            catch
                FRHeatrast{ar,iv} = iFRHRsmzSorted;
            end
        end
    end
end

%% PLOT=====================================================================
if plotfigs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%% SWR %%%%%%
    %% su raster
    if plotEventTrigSU
        if strcmp(eventType, 'swr')
            figname = 'wtrackSWRSU';
        elseif strcmp(eventType, 'lick')
            figname = 'wtrackLickSU';
        end
        Pp = load_plotting_params({'defaults', figname});
        for a = 1:length(F) % per animal
            animal = F(a).animal{3};
            fulltime = F(a).output{1}(1).time';
            sT = knnsearch(fulltime, -Pp.win(1));
            eT = knnsearch(fulltime, Pp.win(2));
            time = fulltime(sT:eT);
            bintime = [fulltime(1):Pp.bin:fulltime(end)]';
            sB = knnsearch(bintime, -Pp.win(1));
            eB = knnsearch(bintime, Pp.win(2));
            for c = 1:length(F(a).output{1}) % per su
                ifig = init_plot(showfigs, Pp.position);
                day = F(a).output{1}(c).index(1);
                eps = F(a).output{1}(c).index(4:end);
                nt = F(a).output{1}(c).index(2);
                clust = F(a).output{1}(c).index(3);
                carea = F(a).output{1}(c).area;
                csubarea = F(a).output{1}(c).subarea;
                fprintf('%s %d %d %d\n', animal, day, nt, clust);
                
                %% spike raster
                sf = subaxis(3, 1, [1 2], Pp.posparams{:});
                sf.Tag = 'raster';
                [xx,yy] = find(F(a).output{1}(c).psth(:,sT:eT)');
                spikeTwin = (xx*1e-3)-Pp.win(1);
                h = scatter(spikeTwin, yy, Pp.spikeSz, '.k', 'markeredgealpha', ...
                    Pp.spikeAlpha);
                axis tight
                xlim([-Pp.win(1) Pp.win(2)])
                xticks([]);
                ylabel('Event #');
                line([0 0], ylim, 'linestyle', '-', 'color', [.5 .5 1 .5], 'linewidth', 2)
                
                %% PSTH
                sf = subaxis(3,1,3,Pp.posparams{:});
                sf.Tag = 'psth';
                [xx,yy] = find(F(a).output{1}(c).psth');
                spikeTwin = (xx*1e-3)-abs(fulltime(1));
                h = histc(spikeTwin, bintime);
                hs = smoothdata(h, 1,'loess',20);
                area(bintime(sB:eB), hs(sB:eB), 'facecolor', 'k')
                axis tight
                xlabel('time s')
                ylabel('count')
                line([0 0], ylim, 'linestyle', '-', 'color', [.5 .5 1 .5], 'linewidth', 2)
                
                %%
                allAxesInFigure = findall(gcf,'type','axes');
                linkaxes(allAxesInFigure, 'x');
                stit = sprintf('%s %s %d %d %d %s %s', figname, animal, day, ...
                    nt, clust, carea, csubarea);
                setSuperAxTitle(stit);
                if pausefigs
                    pause
                end
                if savefigs
                    strsave = save_figure([pconf.andef{4} '/' figname '/' animal], stit);
                end
            end
        end
    end
    
%     %% heatraster per animal
%     if plotSWRTrigHeatrastPerAni
%         figname = 'SWRTrigHeatraster';
%         Pp = load_plotting_params({'defaults', figname}); % load params
%         for a = 1:length(modF) % per animal
%             animal = modF(a).animal{3};
%             for ar = 1:length(Fp.areas) % per area
%                 numESet = length(modF(a).dmatIdx);
%                 % find cells in this area
%                 areaIdx = strcmp({modF(a).output{1}.area}', Fp.areas{ar}{1});
%                 subareaIdx = ~cellfun(@isempty, strfind({modF(a).output{1}.subarea}', ...
%                     Fp.areas{ar}{2}));
%                 iareaIdx = find(all([areaIdx subareaIdx],2));
%                 iareaIdx = iareaIdx(~arrayfun(@(x) isempty(modF(a).output{1}(x).mPctChange), ...
%                     iareaIdx,'un',1));
%                 if isempty(iareaIdx)
%                     continue
%                 end
%                 ifig = init_plot(showfigs, Pp.position); % init fig per area
%                 for iv = 1:numESet % per eventSet in dmat
%                     setID = modF(a).dmatIdx{iv};
%                     sf = subaxis(1, numESet, iv, Pp.posparams{:});
%                     sf.Tag = 'heatraster';
%                     % make firing rate heatraster from all the clusters in this area
%                     iFRHR = cell2mat(arrayfun(@(x) modF(a).output{1}(x).evMean{iv}, ...
%                         iareaIdx,'un',0));
%                     iFRHRzsm = smoothdata(zscore(iFRHR,[],2),2,'loess', 100);
%                     pctChange = arrayfun(@(x) modF(a).output{1}(x).mPctChange{iv}, ...
%                         iareaIdx,'un',1);
%                     [~, srtIdx] = sort(pctChange, 1, 'descend');
%                     iFRHRsmzSorted =  iFRHRzsm(srtIdx,:);                   
%                     imagesc(modF(a).output{1}(iareaIdx(1)).time, 1:length(iareaIdx), ...
%                         iFRHRsmzSorted);
%                    
%                     caxis(sf,'auto')
%                     cm = colormap(1-brewermap(100, 'spectral'));
%                     colormap(cm);
%                     h = colorbar;
%                     ylabel(h, 'zfr')
%                     line([0 0], ylim, 'Color', 'k')
%                     xlabel('time s')
%                 end
%                 stit = sprintf('%s %s %s %s %s', figname, animal, Fp.env, setID, strjoin(Fp.areas{:,ar}));
%                 setSuperAxTitle(stit);
%                 if pausefigs
%                     pause
%                 end
%                 if savefigs
%                     strsave = save_figure([pconf.andef{4} '/' figname '/' animal], stit);
%                     F(a).figs.spikePhaseCumPolar{ar} = strsave;
%                 end
%             end
%         end
%     end
    
    %% heatraster all animals
    % need to trim to bin
    if plotEventTrigHeatrastAllAni
        if strcmp(eventType, 'swr')
            figname = 'wtrackSWRHeatRastAllAn';
        elseif strcmp(eventType, 'lick')
            figname = 'wtrackLickHeatRastAllAn';
        end
        Pp = load_plotting_params({'defaults', figname}); % load params
        time = modF(1).output{1}(1).time';
        s = knnsearch(time, -Pp.win(1));
        e = knnsearch(time, Pp.win(2));
        time = time(s:e);
        for ar = 1:size(FRHeatrast,1)
            ifig = init_plot(showfigs, Pp.position); % init fig per area
            for iv = 1:size(FRHeatrast,2)
                sf = subaxis(1, size(FRHeatrast,2), iv, Pp.posparams{:});
                sf.Tag = 'heatraster';
                imagesc(time, 1:size(FRHeatrast{ar,iv},1), FRHeatrast{ar,iv}(:,s:e));
%                 colormap(magma);
%                 cm = colormap(1-brewermap(100, 'spectral'));
                colormap(viridis);
                h = colorbar;
                caxis(sf,[-4 4])
                ylabel(h, 'zscore firing rate')
                line([0 0], ylim, 'Color', 'k', 'linestyle', '--')
                title(modF(1).dmatIdx{iv})
                xlabel('time s')
                ylabel('Unit #')
            end
            stit = sprintf('%s %s %s', figname, strjoin(Fp.areas{:,ar}));
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure([pconf.andef{4} '/' figname], stit);
            end
        end
    end
    
    %% cdf and stingray
    if plotEventTrigModCDF
        % for each area, plot CDF of eventSets vs shuf
        if strcmp(eventType, 'swr')
            figname = 'wtrackSWRModCDF';
        elseif strcmp(eventType, 'lick')
            figname = 'wtrackLickModCDF';
        end
        Pp = load_plotting_params({'defaults', figname}); % load params
        for ar = 1:size(iMZChSh,1)
            for iv = 1:size(iMPctCh,2)
                ifig = init_plot(showfigs, Pp.position); % init fig per area
                %% cdf
                sf = subaxis(1, 1, 1, Pp.posparams{:});
                sf.Tag = 'ecdf';
                numbins = length(iMPctChSh{ar,iv});
                
                [shufh, shufb] = histcounts(iMPctChSh{ar,iv}, numbins,'Normalization','cdf');
                plot(shufb(1:end-1)+diff(shufb(1:2)) / 2, shufh, 'k')
                hold on;
                [h, b] = histcounts(iMPctCh{ar,iv}, numbins, 'Normalization', 'cdf');
                plot(b(1:end-1) + diff(b(1:2)) / 2, h, 'color', 'b')
%                 set(gca, 'XScale', 'log')
                xlabel('% Change From Baseline')
                ylabel('% units');
                ax = gca;
%                 ax.YDir = 'reverse';
                axis tight
                
                sort_iMPctChSh = sort(iMPctChSh{ar,iv}, 'ascend');
                sigl = sort_iMPctChSh(round(length(iMPctChSh{ar,iv})*.95));
                line([sigl sigl], ylim, 'linestyle', '--', 'color', [.5 .5 .5 .8], 'linewidth', 1)
                hold off

                %% stingray
                
%                 sf = subaxis(1, 2, 2, Pp.posparams{:});
%                 sf.Tag = 'stingray';
%                 
%                 grps = [zeros(numel(iMPctChSh{ar,iv}),1); ones(numel(iMPctCh{ar,iv}),1)];
%                 violin({iMPctChSh{ar,iv}, iMPctCh{ar,iv}},...
%                     'facecolor',[.5 .5 .5; .6 .6 1;],'edgecolor','none');
%                 
%                 legend off
%                 hold on
%                 b = boxplot(sf, [iMPctChSh{ar,iv}; iMPctCh{ar,iv}],grps, ...
%                     'PlotStyle', 'compact', 'Symbol', '.','Color', 'k');
%                     
%                 set(gca,'XTickLabel',{' '})
%                 
%                 xticks([1 2])
%                 xticklabels({'shuffle', 'ca1 response'})
%                 set(gca, 'FontSize', 10)
%                 legend off
%                 hold on;
%                 [p,h,stats] = ranksum(iMPctChSh{ar,iv}, iMPctCh{ar,iv});
%                 xlabel(sprintf('ranksum p%.05f', p),'fontname','arial','fontsize', 10);
%                 ylabel('% from baseline')
%                 hold off;
                
                stit = sprintf('%s %s %s %s', figname, Fp.env, strjoin(Fp.areas{:,ar}), modF(1).dmatIdx{iv});
                setSuperAxTitle(stit);
                if pausefigs
                    pause
                end
                if savefigs
                    strsave = save_figure([pconf.andef{4} '/' figname], stit);
                end
            end
        end
    end
end
%     %% %%%%%% Lick %%%%%%
%     if plotLickTrigSU
%        
%         
%     end
%     %% Plot all lickburst licks. (all an, events) (per area, celltype)
%     if plotLickTrigHeatRaster
%         % time =
%         % mod_Area1_type1_lickTrigSpiking =
%         % mod_Area1_type1_lickTrigSpiking_sh =
%         
%         plot(time, Area1_type1_lickTrigSpiking)
%     end
%     %% sig plot, testing
%     if plotLickTrigModCDF
%         % cdf plot vs shuf
%         [h, e] = histcounts(mod_Area1_type1_lickTrigSpiking, 200, 'Normalization', 'cdf');
%         [hsh, e] = histcounts(mod_Area1_type1_lickTrigSpiking_sh, 200, 'Normalization', 'cdf');
%         
%         plot(e, h);
%         hold on;
%         
%         if pctSU_sigmod_CA1_PN_lickTrigSpiking > .05 %?
%             fprintf('result: area %s is sig mod', iarea)
%         end
%     end
% end

%% %%%%%% ILI Phase %%%%%%
%% all lickburst ILI. (all an, events) (per area, celltype)

                % get signal with uniform prior distribution.
                % compute circular ranks of the phase distribution
%                 y = 
%                 % ecdf
%                 [f,x] = ecdf(y);
%                 % inverse ecdf
%                 unfPhDist = 1/f; 
                
                
                
                % for phasemod: sigmod lines
%                 sigP = [.05]; % .01 .001];
%                 sigL = log(-log(sigP));
%                 line([sigL sigL], ylim, 'linestyle', '--', 'color', 'k');
                
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% LFP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%














