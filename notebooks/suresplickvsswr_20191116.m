%{
NB 20191116
finalize event trig results
.. continued eventtrig_20191117

eventType: lick, ca1rippleskons
eventSet: all, nonlickburst.all, lickburst.all, lickburst.wet, lickburst.dry
data: spikes, eeg
area: ca1, mec
env: wtrack

Figure 2 + supp. characterize event triggered neural activity 
NBGoals: 
1. event trig spiking: What fraction SU (Type,Area) are Modulated
  ++ dfa_eventTrigSpiking (barn.rat) -> calcMod -> 
    - per subject, per eventSet, per area, per su
        -> plotraster?? -> isSig??
    - per subject, per eventSet, per area, per suType
        -> plotHeatRaster?? -> plotCDFsig??
    - all subject, per eventSet, per area, per suType
        -> plotHeatRaster?? -> plotCDFsig??


2. event trig lfp TxF power, itpc: Any Area Sig zmask Mod for power or phase?
 ++ dfa_eventTrigLFP (forest.bear) > stack_riptriglfp > computeAnalyticSignal > getPower >
        - per subject, per eventSet, per area, per nt
                -> licktrigLFPspect_20191103:: plot_expvarCatMeanPwr
        - per subject, per eventSet, per area
                -> licktrigLFPspect_20191103:: combineArea > plot_ByArea
        - all subject, per eventSet, per area
                -> licktrigLFPspect_20191103:: combineArea > ??combineAnimals?? > ??plot??

3. su spike ILI phase clustering (and polar spectrogram?)
- dfa_lickXCorrSpikes (barn.pig) -> plotClustILPC -> plotILPTHeatRaster -> plotCumPolar
    ///alternative strategy transform time in event trig stack to ILI-phase

Notes:
Design Matrix for eventSet. given (an,day)eventtime, return DM
- [make_expvarCat:make_expvarCatMeanPwr]

- *** valid_ntrodes filter is not currently being respected for spiking
%}
pconf = paramconfig;
eventTrigLFP = 0;
eventTrigSpiking = 1;

% run FF
create_filter = 1;
run_ff = 1;
load_ffdata = 0;

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
plot_expvarCatMeanPwr = 0;
plot_ByArea = 0;
pausefigs = 1;
savefigs = 0;

%%
Fp = [];
Fp.animals = {'JZ4'}; %, 'JZ1', 'JZ4'};
eventType = 'swr';
areas = {{'ca1', 'd'}, {'mec', 'deep'}, {'mec', 'supf'}};

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

%% make swr-trig lfp tensor
if stack_LFP
    lfpstack = stack_riptriglfp(F, Fp);
end
if load_LFPstack
    lfpstack = load_data(Fp.paths.resultsDirectory, ['tensor_' Fp.Label], Fp.animals);
end
%% make rawpwr (all trials) [ntrode time rip freq]
if make_rawpwr
    [rawpwr, ~] = computeAnalyticSignal(lfpstack, 'waveSet', Fp.waveSet, 'saveOutput',1, ...
        'lfptype', Fp.uselfptype, 'env', Fp.env, 'eventType', Fp.eventType); % uses parfor
end
if load_rawpwr
    rawpwr = load_data(sprintf('%s/analyticSignal/', pconf.andef{2}), ...
        sprintf('LFPpower_%s_%s_%s_%s', wp.waveSet, Fp.uselfptype, Fp.env, ...
        Fp.eventType), Fp.animals);
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

%% get mean power per condition
if make_expvarCatMeanPwr % :expvarCat @mean /ntTF $time
    evMPwr = getPower(expvarCat, rawpwr, Fp, 'run_perm', 0, 'eventType', Fp.eventType);
end
if load_expvarCatMeanPwr
    outdir = 'expvarCatMeanPwr';
    outpath = [pconf.andef{2},outdir,'/'];
    evMPwr = load_data(outpath, [outdir,'_', Fp.env '_' Fp.eventType], Fp.animals);
end

%% combine perArea perCondition
if combineArea
    for ani = 1:length(evMPwr) % for each animal
        animal = evMPwr(ani).animal;
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
%         ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'') && 'ref'');
%         ntrodes = unique(ntrodes(:,3));
        ntrodes = evMPwr(ani).ntrode;
        for ia = 1:length(areas)
            ntsInArea = evaluatefilter(ntinfo,...
                sprintf('isequal($area,''%s'') && isequal($subarea,''%s'')', areas{ia}{1}, ...
                areas{ia}{2}));
            ntsInArea = unique(ntsInArea(:,3));
            ntsAIdx = knnsearch(ntrodes, ntsInArea);
            for iv = 1:length(evMPwr(ani).expvars)
                areaData = evMPwr(ani).data{iv}.pwr_mean_db(ntsAIdx,:,:);
                evMPwrArea{ia}{iv} = squeeze(nanmean(areaData,1))';
            end
        end
    end
end

%% need to combine animals, per celltype