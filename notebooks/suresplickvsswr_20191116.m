%{
NB 20191116
finalize event trig results

eventType: lick, ca1rippleskons
eventSet: all, nonlickburst.all, lickburst.all, lickburst.wet, lickburst.dry
data: spikes, eeg
area: ca1, mec
env: wtrack

NBGoals: 
1. event trig spiking
    - per subject, per eventSet, per area, per su
        {}
        + dfa_eventTrigSpiking (barn.rat) -> calcMod -> 
    - per subject, per eventSet, per area, per suType
        {}
        + dfa_eventTrigSpiking (barn.rat) -> calcMod -> plotHeatRaster
    - all subject, per eventSet, per area, per suType
        {}
        + dfa_eventTrigSpiking (barn.rat) -> calcMod -> plotHeatRaster

2. pct su (PN,FS) sig modulated, test
    - per subject, per eventSet, per area
    - all subject, per eventSet, per area

3. event trig lfp TxF power, itpc w sig zmask
 ++ dfa_eventTrigLFP (forest.bear) > stack_riptriglfp > computeAnalyticSignal > getPower >
        - per subject, per eventSet, per area, per nt
                -> licktrigLFPspect_20191103:: plot_expvarCatMeanPwr
        - per subject, per eventSet, per area
                -> licktrigLFPspect_20191103:: combineArea > plot_ByArea
        - all subject, per eventSet, per area
                -> licktrigLFPspect_20191103:: combineArea > ??combineAnimals?? > ??plot??

Notes:

- dfa_eventTrigSpiking (barn.rat) -> calcMod -> plotHeatRaster
- dfa_eventTrigLFP (forest.bear) -> stack_LFP -> make_rawpwr -> [make_expvarCat:make_expvarCatMeanPwr]

= transform event trig spike stack to ILI-phase:
- to do? this would be an alternative strategy to using dfa_lickXCorrSpikes

su spike ILI phase clustering
- dfa_lickXCorrSpikes (barn.pig) -> plotClustILPC -> plotILPTHeatRaster -> plotCumPolar

%}
pconf = paramconfig;
% get data
create_filter = 0;
run_ff = 0;
load_ffdata = 0;
% stack data
stack_LFP = 0;
load_LFPstack = 0;
% get power
make_rawpwr = 1;
load_rawpwr = 0;
% create condition design mat
make_expvarCat = 0;
load_expvarCat = 0;
% compute per condition
make_expvarCatMeanPwr = 0;
load_expvarCatMeanPwr = 0;
% combine per area
combineArea = 0;
% plot
plot_expvarCatMeanPwr = 0;
plot_ByArea = 0;
pausefigs = 1;
savefigs = 0;

%% 
Fp.animals = {'JZ4'};
Fp.filtfunction = 'dfa_eventTrigLFP';
Fp.Label = 'wtrackLickTrigLFP';
% {'sleepwtrackdays', 'excludeNoise'};
Fp.params = {'wtrackdays', 'valid_ntrodes', 'lickbouts', 'excludePriorFirstWell', ...
    'excludeAfterLastWell', 'referenced', '4-350Hz',  Fp.Label, Fp.filtfunction};
Fp = load_filter_params(Fp);
wp = getWaveParams(Fp.waveSet);
rs = ''; % ripstate
areas = {{'ca1', 'd'}, {'mec', 'deep'}, {'mec', 'supf'}};

%%
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', ...
        Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
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