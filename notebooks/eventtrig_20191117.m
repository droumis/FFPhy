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
        -> stack_riptriglfp (cactus:: ) ->
        -> computeAnalyticSignal (mushroom:: ) ->
        -> makeEventSet (beer:: ) ->
        -> get_power (cactus:: eventSetPowerTrace) ->

% holy shit.. jz1 licktrigspiking is like 10 times the size on ram than it
is on disk... i need to curate the licks to a subset

also need to filter licks and ILI to those within lickburst range
%}

pconf = paramconfig;
eventTrigLFP = 0; % PIPE:forest.bear.cactus.mushroom.beer.leaf == eventSet mean Spect
eventTrigSpiking = 1; % PIPE:barn.rat.beer.wheelbarrow == eventSet SU mod
eventType = 'lick'; %lick swr

% run FF
create_filter = 0;
run_ff = 0;
load_ffdata = 0;

%% LFP to Analytic Signal Stack
% stack data
stack_LFP = 0; % Comely cactus
load_LFPstack = 0;
% get power, phase of all
make_rawpwr = 0; % Mendacious Mushroom
load_rawpwr = 0;

%% DESIGN MAT MAKER.. WORKS WITH SPIKE AND LFP
make_expvarCat = 0; % Baleful beer
load_expvarCat = 0;

%% LFP Mean Power Per ExpVar
make_expvarCatMeanPwr = 0; % Lachrymose leaf
load_expvarCatMeanPwr = 0;
% combine per area
combineArea = 0;

%% LFP ITPC Per ExpVar

%% Spike Event Time Modulation
calcSUMod = 0; % Wheedling wheelbarrow
loadSUMod = 0;

%% Spike Event PhaseModulation
calcSUphasemod = 0;
loadSUPhaseMod = 0;
%Gather across animals, per area, condition
gatherTimeModResults = 0;
gatherPhaseModResults = 1;

%% SWR x Event Phasic Modulation

%% plot
plotfigs = 0;
showfigs = 1;
pausefigs = 1;
savefigs = 1;

% LFP
plotLFPPerAreaAllAn = 0;

% spike
plotEventSU = 0;         % per eventSet, per SU
plotEventHeatRast = 0; % per eventSet, per area, per animal and all animals
plotEventHeatRastAllAni = 0; % requires gatherResults
plotEventModCDF = 0;     % requires gatherResults// per eventSet, per area, per animal and all animals
% spike phasemod lick
plotSpikePhaseModHeatRast = 0;
%%
Fp = [];
Fp.animals = {'D10', 'JZ1', 'JZ4'};
Fp.areas = {{'ca1', 'd'}, {'mec', 'deep'}, {'mec', 'supf'}};

if eventTrigLFP
    Fp.filtfunction = 'dfa_eventTrigLFP'; % Bellicose Bear
    if strcmp(eventType, 'lick')
        expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
        Fp.Label = 'wtrackLickTrigLFP';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
        'excludeAfterLastWell', 'referenced', '4-350Hz',  Fp.Label, Fp.filtfunction};
    elseif strcmp(eventType, 'swr')
        expvars = {'all', 'lickbouts', 'nolickbouts'};
        Fp.Label = 'wtrackSWRTrigLFP';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
        'excludeAfterLastWell', 'referenced', '4-350Hz',  'ripples', ...
        Fp.Label, Fp.filtfunction};
    end
    
    
elseif eventTrigSpiking
    Fp.filtfunction = 'dfa_eventTrigSpiking'; % Redolent Rat
    if strcmp(eventType, 'lick')
        expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
        Fp.Label = 'wtrackLickTrigSpiking';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
            'excludeAfterLastWell', 'nonMU_cells', Fp.Label, Fp.filtfunction};
    elseif strcmp(eventType, 'swr') %'excludeNoise', 
        expvars = {'all', 'lickbouts', 'nolickbouts'};
        Fp.Label = 'wtrackSWRTrigSpiking';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
            'excludeAfterLastWell', 'nonMU_cells', 'ripples', ...
            Fp.Label, Fp.filtfunction}; % 'excludeNoise',
    end
end
Fp = load_filter_params(Fp);
if eventTrigLFP
    wp = getWaveParams(Fp.waveSet);
end
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

% %% phasic spike mod
% if phaseSpikeMod
%     for a = 1:length(F)
%         animal = F(a).animal{3};
%         for c = 1:length(F(a).output{1})
%             eT = F(a).output{1}(c).eventTimes;
%             ILI = diff(eT);
% %             o(a).output(e).ILI = deT;
%              bin = diff(F(a).output{1}(c).time(1:2)); % should be ~1 ms
%              cIdx = knnsearch(F(a).output{1}(c).time', 0);
%              deTIdx = knnsearch(F(a).output{1}(c).time', ILI);
% %              o(a).output(e).ILIidx = deTIdx;
%             pinceLick = [];
%             for e = 1:length(ILI)
%                 spIli = F(a).output{1}(c).psth(e,cIdx:deTIdx(e));
%                 if any(spIli)
%                     spkOffset = find(spIli)*bin; % idx distance from center (event), scaled to time
%                     pinceLick{e,1} = [spkOffset ./ ILI(e)]';
%                 end
%             end
%             spikePctSinceLick = cell2mat(pinceLick);
%             if ~isempty(spikePctSinceLick)
%                 spikeILIphase = spikePctSinceLick*2*pi;
%                 meanvec = mean(exp(1i*spikeILIphase));
%                 meanMRVmag = abs(meanvec);
%                 vecang = angle(meanvec);
%                 [~, z] = circ_rtest(spikeILIphase);
%                 phasemod = log(z);
%                 F(a).output{1}(c).spikePctSinceLick = spikePctSinceLick;
%                 F(a).output{1}(c).spikeLickPhase = spikeILIphase;
%                 F(a).output{1}(c).meanMRVmag = meanMRVmag;
%                 F(a).output{1}(c).vecang = vecang;
%                 F(a).output{1}(c).phasemod = phasemod;
%             else
%                 F(a).output{1}(c).spikePctSinceLick = [];
%                 F(a).output{1}(c).spikeLickPhase = [];
%                 F(a).output{1}(c).meanMRVmag = [];
%                 F(a).output{1}(c).vecang = [];
%                 F(a).output{1}(c).phasemod = [];
%             end
%         end
%     end
% end

%% vectorize swr-trig lfp
if stack_LFP
    lfpstack = stack_riptriglfp(F, Fp);
end
if load_LFPstack
    lfpstack = load_data(Fp.paths.resultsDirectory, ['tensor_' Fp.Label], Fp.animals);
end

%% make rawpwr (all trials) [ntrode time rip freq]
if make_rawpwr
    % i need to fix this so that it doesn't crash anything but a super computer
    % make it memory aware, and calculate how many workers i can run at
    % once. 
    [rawpwr, ~] = computeAnalyticSignal(lfpstack, 'waveSet', Fp.waveSet, 'saveOutput',1, ...
        'lfptype', Fp.uselfptype, 'env', Fp.env, 'eventType', Fp.eventType); % uses parfor
end
if load_rawpwr
    rawpwr = load_data(sprintf('%s/analyticSignal/', pconf.andef{2}), ...
        sprintf('LFPpower_%s_%s_%s_%s', wp.waveSet, Fp.uselfptype, Fp.env, ...
        Fp.eventType), Fp.animals);
end

%% FOR SPIKES make design mat to slice trials
if make_expvarCat
    if eventTrigSpiking
        data = [];
        for a = 1:length(F)
            data(a).animal = F(a).animal;
            idx = cell2mat({F(a).output{1}.index}');
            [~, dayUnqIdx] = unique(idx(:,1));
            days = idx(dayUnqIdx,[1]);
            eps = idx(dayUnqIdx,[4:5]);
            numEvPerDay = cell2mat({F(a).output{1}(dayUnqIdx).numEventsPerEp}');
            data(a).day = [];
            data(a).epoch = [];
            for d = 1:size(eps,1)
                for e = 1:size(eps,2)
                    data(a).day = [data(a).day; repmat(days(d), numEvPerDay(d,e), 1)];
                    data(a).epoch = [data(a).epoch; repmat(eps(d,e), numEvPerDay(d,e),1)];
                end
            end
            data(a).evStart = cell2mat({F(a).output{1}(dayUnqIdx).eventTimes}');
        end
    else
        data = lfpstack;
    end
    dmat = makeExpvarCatDesignMat(data, expvars, 'eventType', Fp.eventType);
end

if load_expvarCat
    outdir = 'expvarCat';
    outpath = [pconf.andef{2},outdir,'/'];
    dmat = load_data(outpath, [outdir,'_',Fp.env,'_',Fp.eventType], Fp.animals);
end

%% LFP POWER per condition
if make_expvarCatMeanPwr % :expvarCat @mean /ntTF $time
    evMPwr = getPower(dmat, rawpwr, Fp, 'run_perm', 0, 'eventType', Fp.eventType);
end
if load_expvarCatMeanPwr
    outdir = 'expvarCatMeanPwr';
    outpath = [pconf.andef{2},outdir,'/'];
    evMPwr = load_data(outpath, [outdir,'_', Fp.env '_' Fp.eventType], Fp.animals);
end

%% Spike Phase mod
if calcSUphasemod
    pmodF = calcPhaseMod(F, dmat);
    save_data(pmodF, 'results', [Fp.Label '_phasemod']);
end
if loadSUPhaseMod
    pmodF = load_data('results', [Fp.Label '_phasemod'], Fp.animals);
end

%% SPIKE Gather all Phasemod su per area, eventSet
if gatherPhaseModResults
    areaPhasemod = {};
    allAnPmodF = cell2mat(arrayfun(@(x) pmodF(x).output{1}', 1:length(pmodF), 'un', 0)');
    numESet = length(pmodF(1).dmatIdx);
%     iMPctCh = [];
%     iMPctChSh = [];
%     PhModSorted = {};
    for ar = 1:length(Fp.areas) % per area
        % find cells in this area
        areaIdx = strcmp({allAnPmodF.area}', Fp.areas{ar}{1});
        subareaIdx = ~cellfun(@isempty, strfind({allAnPmodF.subarea}', Fp.areas{ar}{2}), 'un', 1);
        iareaIdx = find(all([areaIdx subareaIdx],2));
        iareaIdx = iareaIdx(arrayfun(@(x) ~isempty(allAnPmodF(x).phasemod), ...
            iareaIdx,'un',1));
        
        for iv = 1:numESet % per eventSet    
            gud = [];
            dumpy = [];
            phasemod = [];
            mPctChangeSh = [];
            for i = 1:length(iareaIdx)
                try
                    m = allAnPmodF(iareaIdx(i)).phasemod{iv};
%                     mSh = allAnPmodF(iareaIdx(i)).mPctChangeSh{iv};
                    if ~isempty(m)
                        phasemod = [phasemod; m];
                        mPctChangeSh = [mPctChangeSh; mSh];
                        gud = [gud; i];
                    end
                catch
                    continue
                end
            end
            areaPhasemod{ar,iv} = phasemod; %cell2mat(arrayfun(@(x) allmodF(x).mPctChange{iv}, iareaIdx,'un',0));
%             iMPctChSh{ar,iv} = mPctChangeSh; %cell2mat(arrayfun(@(x) allmodF(x).mPctChangeSh{iv}, iareaIdx,'un',0));
            
            % Firing Rate HeatRaster concat, normalize, smooth, sort
%             iFRHRz = cell2mat(arrayfun(@(x) allAnPmodF(x).evMeanZ{iv}, iareaIdx(gud),'un',0));
%             iFRHRzsm = smoothdata(iFRHRz, 2, 'loess', 10);
% %             mPctChange = cell2mat(arrayfun(@(x) allmodF(x).mPctChange{iv}, iareaIdx,'un',0));
%             [~, srtIdx] = sort(mPctChange, 1, 'descend');
%             iFRHRsmzSorted =  iFRHRzsm(srtIdx,:);
%             try
%                 FRHeatrast{ar,iv} = [FRHeatrast{ar,iv}; iFRHRsmzSorted];
%             catch
%                 FRHeatrast{ar,iv} = iFRHRsmzSorted;
%             end
        end
    end
end

%% calc su modulation
if calcSUMod % wheelbarrow
    modF = calcSUmod(F, dmat);
    save_data(modF, 'results', Fp.Label);
end
if loadSUMod
    modF = load_data('results', Fp.Label, Fp.animals);
end

%% SPIKE Gather all su per area, eventSet
if gatherTimeModResults
    FRHeatrast = {};
    allmodF = cell2mat(arrayfun(@(x) modF(x).output{1}', 1:length(modF), 'un', 0)');
    numESet = length(modF(1).dmatIdx);
    iMPctCh = [];
    iMPctChSh = [];
    PhModSorted = {};
    for ar = 1:length(Fp.areas) % per area
        % find cells in this area
        areaIdx = strcmp({allmodF.area}', Fp.areas{ar}{1});
        subareaIdx = ~cellfun(@isempty, strfind({allmodF.subarea}', Fp.areas{ar}{2}), 'un', 1);
        iareaIdx = find(all([areaIdx subareaIdx],2));
        iareaIdx = iareaIdx(arrayfun(@(x) ~isempty(allmodF(x).mPctChange), ...
            iareaIdx,'un',1));
        
        for iv = 1:numESet % per eventSet    
            % modPctRank
%             iMPRank{ar,iv} = cell2mat(arrayfun(@(x) allmodF(x).modPctRank{iv}, iareaIdx,'un',0));
            gud = [];
            dumpy = [];
            mPctChange = [];
            mPctChangeSh = [];
            for i = 1:length(iareaIdx)
                try
                    m = allmodF(iareaIdx(i)).mPctChange{iv};
                    mSh = allmodF(iareaIdx(i)).mPctChangeSh{iv};
                    if ~isempty(m)
                        mPctChange = [mPctChange; m];
                        mPctChangeSh = [mPctChangeSh; mSh];
                        gud = [gud; i];
                    end
                catch
                    continue
                end
            end
            iMPctCh{ar,iv} = mPctChange; %cell2mat(arrayfun(@(x) allmodF(x).mPctChange{iv}, iareaIdx,'un',0));
            iMPctChSh{ar,iv} = mPctChangeSh; %cell2mat(arrayfun(@(x) allmodF(x).mPctChangeSh{iv}, iareaIdx,'un',0));
%             iMZCh{ar,iv} = cell2mat(arrayfun(@(x) allmodF(x).mZChange{iv}, iareaIdx,'un',0));
%             iMZChSh{ar,iv} = cell2mat(arrayfun(@(x) allmodF(x).mZChangeSh{iv}, iareaIdx,'un',0));
%             iMZResp{ar,iv} = cell2mat(arrayfun(@(x) allmodF(x).mZResp{iv}, iareaIdx,'un',0));
            
            % Firing Rate HeatRaster concat, normalize, smooth, sort
            iFRHRz = cell2mat(arrayfun(@(x) allmodF(x).evMeanZ{iv}, iareaIdx(gud),'un',0));
            iFRHRzsm = smoothdata(iFRHRz, 2, 'loess', 10);
%             mPctChange = cell2mat(arrayfun(@(x) allmodF(x).mPctChange{iv}, iareaIdx,'un',0));
            [~, srtIdx] = sort(mPctChange, 1, 'descend');
            iFRHRsmzSorted =  iFRHRzsm(srtIdx,:);
            try
                FRHeatrast{ar,iv} = [FRHeatrast{ar,iv}; iFRHRsmzSorted];
            catch
                FRHeatrast{ar,iv} = iFRHRsmzSorted;
            end
        end
    end
end

%% LFP combine perArea perCondition
if combineArea
    for ani = 1:length(evMPwr) % for each animal
        try
            animal = evMPwr(ani).animal{3};
        catch
            animal = evMPwr(ani).animal;
        end
        evMPwrArea(ani).animal = animal;
        evMPwrArea(ani).expvars = evMPwr(ani).expvars;
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
%         ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'') && 'ref'');
%         ntrodes = unique(ntrodes(:,3));
        ntrodes = evMPwr(ani).ntrode;
        for ia = 1:length(Fp.areas)
            ntsInArea = evaluatefilter(ntinfo,...
                sprintf('isequal($area,''%s'') && isequal($subarea,''%s'')', Fp.areas{ia}{1}, ...
                Fp.areas{ia}{2}));
            ntsInArea = unique(ntsInArea(:,3));
            ntsAIdx = knnsearch(ntrodes, ntsInArea);
            for iv = 1:length(evMPwr(ani).expvars)
                areaData = evMPwr(ani).data{iv}.pwr_mean_db(ntsAIdx,:,:);
                evMPwrArea(ani).data{ia,iv} = squeeze(nanmean(areaData,1))';
            end
        end
    end
%     evMPwrAreaAllAn = {evMPwrArea.data};
end

%% PLOT=====================================================================
if plotfigs
    %% plot phasemod heatrast per animal
    if plotSpikePhaseModHeatRast
        figname = 'suLickPhaseModHeatRastWtrack';
        Pp=load_plotting_params({'defaults', figname});
        binphase = linspace(0, 2*pi, Pp.numBins);
        binphC = binphase(1:end-1) + diff(1:2)/2;
        x = [binphC]; % binphC+max(binphC)];
        for a = 1:length(pmodF)
            animal = pmodF(a).animal{3};
            fprintf('%s\n', animal);
            for ar = 1:length(Fp.areas)
                ifig = init_plot(showfigs, Pp.position);
                % find cells in this area
                areaIdx = strcmp({pmodF(a).output{1}.area}', Fp.areas{ar}{1});
                subareaIdx = ~cellfun(@isempty, strfind({pmodF(a).output{1}.subarea}', Fp.areas{ar}{2}));
                iareaIdx = find(all([areaIdx subareaIdx],2));
                for iv = 1:length(pmodF(a).dmatIdx)
                % make heatraster from all the clusters in this area,
                % condition
                spikeLickPhaseHist = cell2mat(arrayfun(@(x) ...
                    histcounts(pmodF(a).output{1}(x).spikeLickPhase{iv}, binphase, 'Normalization', ...
                    'probability'), iareaIdx, 'un', 0));
%                 spikeLickPhaseHistShift = circshift(spikeLickPhaseHist, round(size(spikeLickPhaseHist,2)/2), 2);
                [~, imx] = max(spikeLickPhaseHist,[],2);
                %             [~, modsortIdx] = sort(max(spikeLickPhase,[],2), 1, 'descend');
                [~, modsortIdx] = sort(imx, 1, 'descend');
                spikeLickPhaseSort = spikeLickPhaseHist(modsortIdx,:);
                spikeLickPhaseZ = zscore(smoothdata(spikeLickPhaseSort, 2, 'loess', 6), [],2);

                %% all su in area heatraster
                sf = subaxis(1,1,1, Pp.posparams{:});
                sf.Tag = 'heatrast';
                imagesc(x, 1:length(iareaIdx), spikeLickPhaseZ);
                caxis(sf,'auto')
                colormap('jet');
                line([0 0], ylim, 'Color', 'k')
                ylabel(sprintf('%s SU num', strjoin(Fp.areas{:,ar})))
                %                 title(sprintf('%s %s', Fp.areas{ar}{1}, Fp.areas{ar}{2}))
                xlabel('lick phase rad')
                %%
                %             allAxesInFigure = findall(gcf,'type','axes');
                %             linkaxes(allAxesInFigure, 'x');
                stit = sprintf('%s %s %s %s', figname, animal, Fp.env, strjoin(Fp.areas{:,ar}));
                setSuperAxTitle(stit);
                if pausefigs
                    pause;
                end
                if savefigs
                    save_figure([pconf.andef{4} '/' figname '/' animal], stit);
                end
                end
            end
        end
    end
    
    %% plot LFP PWR per area per condition
    if plotLFPPerAreaAllAn
        if strcmp(eventType, 'swr')
            figname = 'wtrackSWRPwrAreaPerAn';
        elseif strcmp(eventType, 'lick')
            figname = 'wtrackLickPwrAreaPerAn';
        end
        Pp=load_plotting_params({'defaults',figname});
        for ani = 1:length(evMPwrArea) % for each animal
            animal = evMPwrArea(ani).animal;
            for ia = 1:length(Fp.areas)
                ifig = init_plot(showfigs, Pp.position);
                for iv = 1:length(evMPwrArea(ani).expvars)
                    sf = subaxis(1,length(Fp.areas), iv);
                    frequency = round(wp.frex);
                    time = wp.win(1):1/(wp.srate/wp.dsamp):wp.win(2);
                    trim = knnsearch(time',Pp.win(1)):knnsearch(time', Pp.win(2));
                    ttime = time(trim);
                    tdata = evMPwrArea(ani).data{ia,iv}(:,trim);
                    contourf(ttime,frequency,tdata,40,'linecolor','none')
                    set(gca,'ydir','normal','yscale','log');
                    
                    colormap(Pp.usecolormap)
                    if iv == 1
                        caxis(sf, 'auto')
                        cax = [-max(abs(caxis)) max(abs(caxis))];
                    end
                    caxis(sf, cax) % normalize to the 'all' condition
                    
                    
%                     if ia == 1
%                         caxis([-4 4]); %cax(2)])
%                     else
%                         caxis([-1.5 1.5]); %cax(2)])
%                     end
                    ytickskip = 2:4:wp.numfrex;
                    set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', Pp.tickFsize)
                    title(sprintf('%s',evMPwr(ani).expvars{iv}))
                    yl = ylim;
                    xl = xlim;
                    line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
                    colorbar
                    hold on
                end
%                 ylabel(sf{1},'frequency Hz')
%                 xlabel(sf{1}, 'time s')
                %%
                stit = sprintf('%s %s %s %s %s %s %s',figname, Fp.areas{ia}{1},Fp.areas{ia}{2}, Fp.eventType, ...
                    animal, Fp.env, Fp.waveSet);
                setSuperAxTitle(stit);
                if pausefigs
                    pause
                end
                if savefigs
                    save_figure([pconf.andef{4} figname], stit);
                end
            end
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%% SWR %%%%%%
    %% su raster
    if plotEventSU
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
                hs = smoothdata(h, 1,'loess', 20);
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
    
    %% heatraster per animal
    if plotEventHeatRast
        if strcmp(eventType, 'swr')
            figname = 'wtrackSWRHeatrast';
        elseif strcmp(eventType, 'lick')
            figname = 'wtrackLickHeatRast';
        end
        Pp = load_plotting_params({'defaults', figname}); % load params
        for a = 1:length(modF) % per animal
            animal = modF(a).animal{3};
            for ar = 1:length(Fp.areas) % per area
                numESet = length(modF(a).dmatIdx);
                % find cells in this area
                areaIdx = strcmp({modF(a).output{1}.area}', Fp.areas{ar}{1});
                subareaIdx = ~cellfun(@isempty, strfind({modF(a).output{1}.subarea}', ...
                    Fp.areas{ar}{2}));
                iareaIdx = find(all([areaIdx subareaIdx],2));
                iareaIdx = iareaIdx(~arrayfun(@(x) isempty(modF(a).output{1}(x).mPctChange), ...
                    iareaIdx,'un',1));
                if isempty(iareaIdx)
                    continue
                end
                ifig = init_plot(showfigs, Pp.position); % init fig per area
                for iv = 1:numESet % per eventSet in dmat
                    setID = modF(a).dmatIdx{iv};
                    sf = subaxis(1, numESet, iv, Pp.posparams{:});
                    sf.Tag = 'heatraster';
                    mPctChange = [];
                    dumpy = [];
                    gud = [];
                    for i = 1:length(iareaIdx)
                        try
                            if ~isempty(modF(a).output{1}(iareaIdx(i)).mPctChange{iv})
                                
                                mPctChange = [mPctChange; modF(a).output{1}(iareaIdx(i)).mPctChange{iv}];
                                gud = [gud; i];
                                
                            else
                                dumpy = [dumpy; i];
                            end
                        catch
                            dumpy = [dumpy; i];
                        end
                    end
                    iareaIdx = iareaIdx(gud);
                    % make firing rate heatraster from all the clusters in this area
                    iFRHRz = cell2mat(arrayfun(@(x) modF(a).output{1}(x).evMeanZ{iv}, ...
                        iareaIdx,'un',0));
                    iFRHRzsm = smoothdata(iFRHRz,2,'loess', 10);
                    
%                     pctChange = cell2mat(arrayfun(@(x) modF(a).output{1}(x).mPctChange{iv}, ...
%                         iareaIdx,'un',0));
                    [~, srtIdx] = sort(mPctChange, 1, 'descend');
                    iFRHRsmzSorted =  iFRHRzsm(srtIdx,:);       
                    time = modF(a).output{1}(iareaIdx(1)).time';
                    s = knnsearch(time, -Pp.win(1));
                    e = knnsearch(time, Pp.win(2));
                    time = time(s:e);
                    imagesc(time, 1:length(iareaIdx), iFRHRsmzSorted(:,s:e));
                   
                    colormap(viridis);
                    h = colorbar;
                    caxis(sf,[-3 3])
                    ylabel(h, 'zscore firing rate')
                    line([0 0], ylim, 'Color', 'k', 'linestyle', '--')
                    title(modF(a).dmatIdx{iv})
                    xlabel('time s')
                    ylabel('Unit #')
                end
                stit = sprintf('%s %s %s %s', figname, animal, strjoin(Fp.areas{:,ar}));
                setSuperAxTitle(stit);
                if pausefigs
                    pause
                end
                if savefigs
                    strsave = save_figure([pconf.andef{4} '/' figname '/' animal], stit);
                    F(a).figs.spikePhaseCumPolar{ar} = strsave;
                end
            end
        end
    end
    
    %% heatraster all animals
    % need to trim to bin
    if plotEventHeatRastAllAni
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
    if plotEventModCDF
        % for each area, plot CDF of eventSets vs shuf
        if strcmp(eventType, 'swr')
            figname = 'wtrackSWRModCDF';
        elseif strcmp(eventType, 'lick')
            figname = 'wtrackLickModCDF';
        end
        Pp = load_plotting_params({'defaults', figname}); % load params
        for ar = 1:size(iMPctChSh,1)
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

%% %%%%%% ILI Phase %%%%%%
%% all lickburst ILI. (all an, events) (per area, celltype)
% pm = [F(a).output{1}.phasemod];
% 

[h, b] = histcounts(areaPhasemod{1,1}, 200, 'Normalization', 'cdf');
plot(b(1:end-1) + diff(b(1:2)) / 2, h)
xlabel('phasemod')
ylabel('% units');
ax = gca;
ax.YDir = 'reverse';
axis tight
hold on

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














