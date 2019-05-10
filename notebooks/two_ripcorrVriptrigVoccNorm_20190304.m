
%{
- fix results saving to make it one per animal
- remake ripcorr, riptrig, behave, ratemap and save per animal
- make jz2, d12 riptrig spiking
- remake riptrig spiking and include all epochs
- compute corrcoef per swr in new analysis function
- load ripcorr, riptrig, behave, ratemap

- combine results per day, subplots tight, perform with conf bounds
- next: compute per rip measure
- next: compute autocorr w score
%}
%% ====RIP CORR================================
% Fp = load_filter_params('mua_calcxcorrmeasures');

%% ====RIP TRIG================================
% Fp = load_filter_params('riptrigspiking_withMU');
% running that again is going to take too long.. instead load the joint data and
% save out by animal
me = animaldef('demetris');
riptrig_filename = 'D10-D13-JZ1-JZ3-JZ4_dfa_riptrigspikings_wtrack.mat';
riptrig_result = load(sprintf('%s/dfa_riptrigspiking/%s', me{2}, ...
    riptrig_filename));
% save_filter_output(riptrig_result.F, paths.filtOutputDirectory, ...
% paths.filenamesave, 'per_animal', 1)

%% ====RATE MAPS====================================
Fp = load_filter_params({'wtrack', 'nonref_ntrodes', ...
    'nonMU_cells','ratemaps'});
Fp.animals = {'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
F = createfilter('animal', Fp.animals, 'epochs', ...
    Fp.epochfilter, 'tetrodes', Fp.tetfilter, 'cells', Fp.cellfilter, ...
    'iterator', Fp.iterator, 'excludetimefilter', Fp.timefilter);
F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes);


%% ====Performance====================================
animals = {'D10'};%, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
for ani = 1:length(animals)
    bstate = load(sprintf(...
        '/data2/demetris/%s/filterframework/%sBehaveState.mat',animals{ani},...
        animals{ani}));
    BehaveState(ani).animal = animals(ani);
    BehaveState(ani).statechanges = bstate.BehaveState.statechanges;
    BehaveState(ani).statespace = bstate.BehaveState.statespace;
%     stspace{ani} = bstate.BehaveState.statespace.allbound;
%     stspace_times{ani} = bstate.BehaveState.statespace.allepsMat(:,3);
end
save_filter_output(BehaveState, ...
    '/data2/demetris/analysis/filter_output/behave_state/', ...
    'behave_state_wtrack', 'per_animal', 1)

%% ====================================
runFilterFramework = 1;
saveFilterOutput = 1;
loadFilterOutput = 0;
%% ---------------- Paths ---------------------------------------------------
paths = make_paths(Fp);
%% ---------------- Run FIlter -----------------------------------------------
if runFilterFramework == 1
    F = runfilter(F);
    for a = 1:length(F)
        F(a).filterparams = Fp;
    end
end
%% ---------------- Save Filter Output ----------------------------------------
if saveFilterOutput == 1;
    save_filter_output(F, paths.filtOutputDirectory, paths.filenamesave, 'per_animal', 1)
end
% 
% loadresults = 1;
% combine_results = 0;
% compute_across_results = 0;
% plotfigs = 0;
% 
% me = animaldef('demetris');
% animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
% 
% if loadresults
%     for an = 1:length(animals);
%         % need to save the results out per animal to make this cleaner.. 
%         % need to save the perform out one combined across days
%         ripcorr_filename = 'D10-D13-JZ1-JZ3-JZ4_mua_calcxcorrmeasuress_wtrack';
%         perform_filename = 
%         riptrig_filename = 'D10-D13-JZ1-JZ3-JZ4_dfa_riptrigspikings_wtrack.mat';
%         ratemap_filename = 
% 
%         riptrig_result{an} = load(sprintf('%s/dfa_riptrigspiking/%s', me{2}, filename));
%         ripcorr_result{an} = load(sprintf('%s/mua_calcxcorrmeasures/%s', me{2}, filename));
%     end
% end
% 
% if combine_results
%     % combine the results for each day
%     % should i make a new combined result or 1 for each result type ?
% end
% 
% if compute_across_results
%     % compute corr coef between rip corr and performance/change
% end
% 
% if plotfigs
%     %% animal loop
%     % for each animal, create a graph for the corr of the rip corr measure Vs ...
%     % performance, per day
%     % could also compute other across-days scores for rip trig, ratemaps,
%     % autocorr, etc
%     
%     %% figure (nt-pair) loop
% end
% 
% 
