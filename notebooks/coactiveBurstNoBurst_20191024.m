

create_filter = 0;
run_ff = 0;
save_ffdata = 0;
load_ffdata = 1;

%% get the lickphase for each swr
% data filter params
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_suCoactiveXcorr';
% condition = 'lickburst';
condition = 'noburst';
% Fp.params = {'savefigs', 'wtrackdays', 'ripples', 'lickbouts', Fp.filtfunction};
Fp.params = {'savefigs', 'wtrackdays', 'ripples', 'nolickbouts', Fp.filtfunction};
% FF
Fp = load_filter_params(Fp);
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'cellpairs', Fp.cellpairfilter, 'iterator',Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff; F = runfilter(F);
    for d = 1:length(F); F(d).datafilter_params = Fp; end
end
if save_ffdata
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', sprintf('_%s_%s', Fp.epochEnvironment, condition));
end
if load_ffdata
    Flb = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', sprintf('_%s_%s', Fp.epochEnvironment, 'lickburst'));
    Fno = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', sprintf('_%s_%s', Fp.epochEnvironment, 'noburst'));
end
%%
% i need to combine epoch results for each cell pair (take mean?)
% plot the cumulative coactive z for each condition. run kstest?
ani = 1;
Flb_struct = cell2mat(Flb(ani).output{1}');
Flb_coZ = [Flb_struct.coactivez]';

Fno_struct = cell2mat(Fno(ani).output{1}');
Fno_coZ = [Fno_struct.coactivez]';
%%
[lickh, lickb] = histcounts(Flb_coZ,1000,'Normalization', 'cdf');
plot(lickb(1:end-1)+diff(lickb) / 2, lickh, 'k')
hold on;
[noh, nob] = histcounts(Fno_coZ,1000,'Normalization', 'cdf');
plot(nob(1:end-1)+diff(nob) / 2, noh, 'b')
hold off
[H, pValue, KSstatistic] = kstest2(Fno_coZ, Flb_coZ);
%%
        
histogram(Flb_coZ, 40, 'Normalization', 'probability')
hold on;
histogram(Fno_coZ, 40, 'Normalization', 'probability', 'DisplayStyle','stairs', 'linewidth', 2)
hold off
legend({'lickburst swr', 'nonburst swr'})
title('coactiveZ ca1')
save_figure([pconf.andef{4} '/coactiveZburstNoburst/'], 'coactiveZburstNoburst_ca1');

% then for each pair get diff of coactiveZ for burst-noburst... and then 
% generate a condition-shuffled distribution of diff coactiveZ.. to
% compare that to?
