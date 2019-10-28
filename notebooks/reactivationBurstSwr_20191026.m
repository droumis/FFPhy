

%{
David wants to see that the swr during lick bursts are meaningful events
that aren't just a consequence of motoric action.. so showing that the
swr's are real by comparing the neural activity during the burst swr to
non-swr times in the lick burst.. and further, that there is some structure
during burst swr's that is different than the activity during nonswr bursts
at the same ILI phase.. 

- compare intra-burst SU PSwrTH to PLickTH (basically just a SU modulation).. 
- compare intra-burst Reactivation PSwrTH to PLickTH

- i think the easier first step is to do the psth relative to time rather
than ILIphase.. 
%}
create_filter = 0;
run_ff = 1;
load_ffdata = 0;
%% FF
Fp.animals = {'D10', 'D12'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_reactivationPSTH';
Fp.params = {'>4cm/s', 'ca1SU', 'wtrackdays', 'exemplar_wepochs', 'excludePriorFirstWell', ...
    Fp.filtfunction};
Fp = load_filter_params(Fp);
%%
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'cells', Fp.cellfilter, 'iterator',Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F), 'un', 1);
    F = runfilter(F);
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.env]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', ['_' Fp.env]);
end



















