
create_filter = 1;
run_ff = 0;
save_ffdata = 0;
load_ffdata = 0;
stack_spikes = 0;

make_licks = 1;

pilotEpoch = 0;
loadEpoch = 0;
run_dfa = 0;

plotfigs = 1;
displayplots = 0;
saveplots = 1;
%% data filter params
Fp.animals = {'JZ3'};
Fp.filtfunction = 'dfa_lickXCorrSpikes';
Fp.add_params = {'savefigs', 'wtrackdays', 'valid_ntrodes', 'nonMU_cells'}; %'exemplar_wepochs'
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
        sprintf('_%s', Fp.epochEnvironment));
end
if load_ffdata
    spikesF = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s', Fp.epochEnvironment));
end

if make_licks
    for ani = 1:length(F)
        animdef = animaldef(F(ani).animal{3});
        DIO = loaddatastruct(animdef{2}, F(ani).animal{3}, 'DIO');
        task = loaddatastruct(animdef{2}, F(ani).animal{3}, 'task');
        get_licks(F(ani).animal{3}, F(ani).epochs{1}, DIO, task)
    end
end
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
        dfa_lickXCorrSpikes([day epoch ntrode cluster], excludeIntervals, 'lick', lick, ...
            'task', task, 'spikes', spikes, 'cellinfo', cellinfo, ...
            'tmax', Fp.tmax, 'bin', Fp.bin, 'plotfigs', plotfigs, 'saveplots', saveplots, ...
            'displayplots', displayplots, 'animal', animal);
    end
end