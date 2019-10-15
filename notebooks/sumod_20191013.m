

%{
previous work in /home/droumis/Src/Matlab/filterframework_dr/notebooks/swrtrigSUmod_20191007.m

all swrs:
    get examples of ca1 swrtrig su rasters, psth
    get ca1 all cells heatraster
    get modscores pop %
lickbout swrs:
    get examples of ca1 swrtrig su rasters, psth
    get ca1 all cells heatraster
    get modscores pop %

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













