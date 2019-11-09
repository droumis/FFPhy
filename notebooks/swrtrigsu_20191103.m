
%{

%}

% get data
pconf = paramconfig;
create_filter = 0;
run_ff = 0;
load_ffdata = 0;
% stack data
combineEpochs = 0;
saveCombinedEpochs = 0;
loadCombinedEpochs = 0;
% get event response
calcmod = 0;
loadmodF = 0;
% plot
plotClusterFigs = 0;
plotPopulationFigs = 0;

pausefigs = 0;
savefigs = 1;

%%
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_eventTrigSpiking';
Fp.Label = 'swrTrigSpiking';
Fp.params = {'ripples', 'exemplar_wepochs', 'valid_ntrodes', 'nonMU_cells', Fp.filtfunction};
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
%     save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
%         'filetail', ['_' Fp.env '_' Fp.eventSourceArea Fp.eventtype]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', ['_' Fp.env '_' Fp.eventSourceArea Fp.eventtype]);
end
%% ---------------- combine cluster epochs ----------------------------
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
    
    modF = calcSUmod(F, filtF, 'respwin', respwin, 'basewin', basewin, 'minNumSwr', ...
        minNumSwr, 'nshuffs', nshuffs, 'shuffms', shuffms, ...
        'filetail', ['_' Fp.epochEnvironment]);
end
% ---------------- load mod---------------------------------------
if loadmodF
    modF = load_data([pconf.andef{3} '/sumod'], 'sumod', Fp.animals,...
        'filetail', ['_' Fp.epochEnvironment]);
end
