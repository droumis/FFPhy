
create_filter = 0;
run_ff = 1;
save_ffdata = 0;
load_ffdata = 0;
%% data filter params
Fp.animals = {'JZ1'}; %, 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_lickswrcorr';
Fp.params = {'pausefigs', 'wtrackdays', Fp.filtfunction}; %'exemplar_wepochs',
%% FF
Fp = load_filter_params(Fp);
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'iterator',Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff; F = runfilter(F);
    for d = 1:length(F); F(d).datafilter_params = Fp; end
end
if save_ffdata
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', sprintf('_%s_%s', Fp.epochEnvironment));
end
if load_ffdata
    lick = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', sprintf('_%s_%s', Fp.epochEnvironment, conditions{1}));
end