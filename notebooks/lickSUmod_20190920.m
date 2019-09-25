    
%{

- for each DENC
- for each lick bout, return a mean vector magnitude 
    - get the phase clustering of the spikes relative to the lick cycle within
        each lick bout
- for each lick bout, randomly scatter the spikes and get a shuffled mean
    vector magnitude

- because of the variability/doublets of the input triggers.. i think for
now the best that could be done is spikes/swr - lfp phase locking during
licking bouts vs outside of licking bouts

- does it make sense to align data to the bout start?
- then do morelet wavelet convolution? it looks like the bouts typically
    last ~ 10 seconds so maybe having -5 s : 15 s... although this is going
    to be kinda weird to show because the animals have high theta phase
    locking while running and often thet trigger the well as soon as they
    reach the well.. so the pre-bout period would have high theta... 
- i guess i'd need to use whatever non-bout immobility periods there are as
the baseline/control condition.. which more likely than not happen after
the active licking..

%}

create_filter = 0;
run_ff = 1;
save_ffdata = 0;
load_ffdata = 0;
stack_spikes = 0;

make_licks = 0;

pilotEpoch = 0;
loadEpoch = 0;
run_dfa = 1;

plotfigs = 1;
displayplots = 0;
saveplots = 1;

%% data filter params
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_lickphaseSUclustering';
Fp.add_params = {'savefigs', 'wtrackdays', 'valid_ntrodes', 'nonMU_cells', ...
    'exemplar_wepochs'};
%% FF
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'iterator',Fp.iterator,'eegtetrodes',Fp.tetfilter,...
        'cells', Fp.cellfilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff; F = runfilter(F); for d = 1:length(F); F(d).datafilter_params = Fp; end; end
if save_ffdata
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s', Fp.epochEnvironment)); end
if load_ffdata
    spikesF = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s', Fp.epochEnvironment)); end

%% pilot one epoch.. bypasses iterator
if pilotEpoch
    day = F(1).epochs{1}(1,1);
    epoch = F(1).epochs{1}(1,2);
    ntrode = F(1).data{1}{1}(end,1);
    cluster = F(1).data{1}{1}(end,2);
    animal = Fp.animals{1};
    animdef = animaldef(animal);
    if loadEpoch
        spikes = loaddatastruct(animdef{2}, animal, 'spikes', day);
        cellinfo = loaddatastruct(animdef{2}, animal, 'cellinfo', day);
        lick = loaddatastruct(animdef{2}, animal, 'lick', day);
        task = loaddatastruct(animdef{2}, animal, 'task', day);
    end
    if run_dfa
        epidx = find(ismember(F.epochs{1}, [day epoch], 'rows'));
        excludeIntervals = F.excludetime{1}{epidx};
        dfa_lickphaseSUclustering([day epoch ntrode cluster], excludeIntervals, 'lick', ...
            lick, 'task', task, 'spikes', spikes, 'cellinfo', cellinfo, 'animal', animal);
%         'tmax', Fp.tmax, 'bin', Fp.bin, 'plotfigs', plotfigs, 'saveplots', saveplots, ...
%             'displayplots', displayplots, 
    end
end