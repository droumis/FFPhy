

%{
- calculate mu corr measures per rip

%}

Fp = load_filter_params({'wtrack','ripples','dfa_perripspikingcorr'});
Fp.animals = {'JZ4'};%, 'JZ4'}; %, , 'JZ1', 'JZ3', 'JZ4'
Fp.days = [1:7];

runFilterFramework = 1;
saveFilterOutput = runFilterFramework;
loadFilterOutput = 0;
load_other_data = 0;

plotstuff = 0;
savefigs = 0;
pausefigs = 0;

%% ---------------- Paths ---------------------------------------------------
paths = make_paths(Fp.filtfunction, Fp.epochEnvironment);
%% ---------------- Run FIlter -----------------------------------------------
if runFilterFramework == 1
    F = createfilter('animal', Fp.animals, 'days', Fp.days,'epochs', ...
        Fp.epochfilter, 'tetrodepairs', Fp.tetpairfilter, 'iterator', ...
        Fp.iterator, 'excludetimefilter', Fp.timefilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes);
    for a = 1:length(F) % save filter detailes along with results
        F(a).datafilter_params = Fp;
    end
    F = runfilter(F);
end
%% ---------------- Save Filter Output ----------------------------------------
if saveFilterOutput == 1;
    save_filter_output(F, paths.filtOutputDirectory, paths.filenamesave)
end
%% ---------------- Load Filter Output ----------------------------------------
if loadFilterOutput == 1;
    F = load_filter_output(paths.filtOutputDirectory, paths.filenamesave);
end