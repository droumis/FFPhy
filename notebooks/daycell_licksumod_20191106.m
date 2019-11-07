
% get data
pconf = paramconfig;
create_filter = 0;
run_ff = 1;
load_ffdata = 0;
% plot
plotfigs = 0;
pausefigs = 1;
savefigs = 0;

%% filter params
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_lickXCorrSpikes';
Fp.Label = 'lickSpikeXcorr';
Fp.params = {'wtrack', '<4cm/s', 'wtrackdays', 'valid_ntrodes', 'nonMU_cells', Fp.filtfunction};
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
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.env '_' Fp.eventType]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', ['_' Fp.env '_' Fp.eventType]);
end