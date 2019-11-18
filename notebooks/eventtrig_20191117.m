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
eventType = 'swr'; %lick

% run FF
create_filter = 0;
run_ff = 0;
load_ffdata = 0;

%% spike
calcMod = 0; % wheelbarrow
loadMod = 0;

%% LFP
% stack data
stack_LFP = 0;
load_LFPstack = 0;
% get power
make_rawpwr = 0;
load_rawpwr = 0;
% create condition design mat
make_expvarCat = 0;
load_expvarCat = 0;
% compute per condition
make_expvarCatMeanPwr = 0;
load_expvarCatMeanPwr = 0;
% combine per area
combineArea = 0;

%% plot
plotfigs = 1;
showfigs = 1;
pausefigs = 1;
savefigs = 1;

% swr
plotSWRTrigSU = 0;         % per eventSet, per SU
plotSWRTrigHeatraster = 1; % per eventSet, per area, per animal and all animals
plotSWRTrigModCDF = 0;     % per eventSet, per area, per animal and all animals

% lick
plotLickTrigSU = 0;         % per eventSet, per SU
plotLickTrigHeatraster = 0; % per eventSet, per area, per animal and all animals
plotLickTrigModCDF = 0;     % per eventSet, per area, per animal and all animals

%%
Fp = [];
Fp.animals = {'D10', 'JZ1', 'JZ4'}; %, 'JZ1', 'JZ4'};
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
            'excludeAfterLastWell', 'nonMU_cells', 'excludeNoise', Fp.Label, Fp.filtfunction};
    elseif strcmp(eventType, 'swr')
        Fp.Label = 'wtrackSWRTrigSpiking';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
            'excludeAfterLastWell', 'nonMU_cells', 'ripples', 'excludeNoise', ...
            Fp.Label, Fp.filtfunction};
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
if calcMod % wheelbarrow
    modF = calcSUmod(F);
    save_data(modF, 'results', Fp.Label);
end
if loadMod
    modF = load_data('results', Fp.Label, Fp.animals);
end

%% PLOT=====================================================================
if plotfigs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %%%% SWR %%%%%%
    %% su raster
    if plotSWRTrigSU
        
    end
    %% heatraster
    if plotSWRTrigHeatraster
        figname = 'SWRTrigHeatraster';
        Pp = load_plotting_params({'defaults', figname}); % load params
        AllAniFRHRsmzSorted = [];
        for a = 1:length(F) % per animal
            animal = F(a).animal{3};
            for ar = 1:length(Fp.areas) % per area
                numESet = length(modF(a).dmatIdx);
                % find cells in this area
                areaIdx = strcmp({F(a).output{1}.area}', Fp.areas{ar}{1});
                subareaIdx = ~cellfun(@isempty, strfind({F(a).output{1}.subarea}', Fp.areas{ar}{2}));
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
                    % make firing rate heatraster from all the clusters in this area
                    iFRHR = cell2mat(arrayfun(@(x) modF(a).output{1}(x).evMean{iv}, ...
                        iareaIdx,'un',0));
                    iFRHRsmz = smoothdata(zscore(iFRHR,[],2),2,'loess', 20);
                    pctChange = arrayfun(@(x) modF(a).output{1}(x).mPctChange{iv}, ...
                        iareaIdx,'un',1);
                    [~, srtIdx] = sort(pctChange, 1, 'descend');
                    iFRHRsmzSorted =  iFRHRsmz(srtIdx,:);
                    AllAniFRHRsmzSorted{a}{ar,iv} = iFRHRsmzSorted;
                    imagesc(modF(a).output{1}(iareaIdx(1)).time, 1:length(iareaIdx), ...
                        iFRHRsmzSorted);
                   
                    caxis(sf,'auto')
%                     cm = colormap(1-brewermap(100, 'spectral'));
                    colormap(magma);
                    h = colorbar;
                    ylabel(h, 'foo')
                    line([0 0], ylim, 'Color', 'k')
                    xlabel('time s')
                end
                stit = sprintf('%s %s %s %s %s', figname, animal, Fp.env, setID, strjoin(Fp.areas{:,ar}));
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
    %% cdf and stingray
    if plotSWRTrigModCDF
        
    end
    
    %% %%%%%% Lick %%%%%%
    %% Plot all lickburst licks. (all an, events) (per area, celltype)
    if plotLickTrigHeatRaster
        % time =
        % mod_Area1_type1_lickTrigSpiking =
        % mod_Area1_type1_lickTrigSpiking_sh =
        
        plot(time, Area1_type1_lickTrigSpiking)
    end
    %% sig plot, testing
    if plotLickTrigModCDF
        % cdf plot vs shuf
        [h, e] = histcounts(mod_Area1_type1_lickTrigSpiking, 200, 'Normalization', 'cdf');
        [hsh, e] = histcounts(mod_Area1_type1_lickTrigSpiking_sh, 200, 'Normalization', 'cdf');
        
        plot(e, h);
        hold on;
        
        if pctSU_sigmod_CA1_PN_lickTrigSpiking > .05 %?
            fprintf('result: area %s is sig mod', iarea)
        end
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% ILI Phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% all lickburst ILI. (all an, events) (per area, celltype)

















