%% Run/Plot rip triggered spiketrain correlation across clusters 

close all
%% ---------------- Dashboard --------------------------

plotfigs = 0;
savefigs = 0;
pausefigs = ~savefigs;
runAntiAlias = 0;

Fp = load_filter_params('riptrigspiking_withMU');
Fp.animals = {'D13'};
Fp.days = [1:7];

runFilterFramework = 1;
loadFilterOutput = 0;

if runFilterFramework || loadFilterOutput
    saveFilterOutput = runFilterFramework;
    combineEpochs = runFilterFramework;
    saveCombinedEpochs = combineEpochs;
    loadCombinedEpochs = loadFilterOutput;

    [F, ppF] = dfs_riptrigspiking(Fp, 'runFilterFramework', runFilterFramework, ...
    'saveFilterOutput', saveFilterOutput, 'loadFilterOutput', loadFilterOutput, ...
    'combineEpochs', combineEpochs, 'saveCombinedEpochs', saveCombinedEpochs, ...
    'loadCombinedEpochs', loadCombinedEpochs);
end

paths = make_paths(Fp);