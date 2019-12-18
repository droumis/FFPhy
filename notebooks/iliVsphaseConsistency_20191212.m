%{
does the swr ILI-phase change along with ILI duration?

- % ILI for each swr vs ILI
%}

pconf = paramconfig;
create_filter = 0;
run_ff = 0;
load_ff = 0;

plotfigs = 0;
showfigs = 0;
pausefigs = 0;
savefigs = 1;
savefigas = {'png', 'eps'};

%% FF Data
Fp = [];
Fp.animals = {'D10', 'JZ1', 'JZ4'};
Fp.areas = {{'ca1', 'd'}, {'mec', 'deep'}, {'mec', 'supf'}};

Fp.filtfunction = 'dfa_lickswrcorr'; % city.alien
Fp.Label = 'wXPTrigSWR';
Fp.params = {'wtrackdays', 'excludePriorFirstWell', 'excludeAfterLastWell', ...
    'ripples>2', Fp.Label, Fp.filtfunction};
Fp = load_filter_params(Fp);

if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'cells', Fp.cellFilter, 'iterator',Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F), 'un', 1);
    F = runfilter(F);
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.env]);
end
if load_ff
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', ['_' Fp.env]);
end

%%
if 1
    swrTSL = {};
    swrPSL = {};
    swrBinILI = {};
    swrBS = {};
    swrLPh = {};
    swrBLN = {};
    swrBCS = {};
    for a = 1:length(F)
        swrBLN{a} = cell2mat({F(a).output{1}.swrBurstLickNum}');
        
%         histogram(swrBLN{2},25)
%         swrBCS{a} = cell2mat(cellfun(@(x) x.burstContainSwr, F(a).output, 'un', 0)');
%         swrLPh{a} = cell2mat(cellfun(@(x) x.swrLickPhase, F(a).output, 'un', 0)');
%         swrBS{a} = cell2mat(cellfun(@(x) x.swrInBurstStart, F(a).output, 'un', 0)');
%         swrPSL{a} = cell2mat(cellfun(@(x) x.swrPctSinceLick, F(a).output, 'un', 0)');
%         swrTSL{a} = cell2mat(cellfun(@(x) x.swrTimeSinceLick, F(a).output, 'un', 0)');
%         swrBinILI{a} = cell2mat(cellfun(@(x) x.swrBinILI, F(a).output, 'un', 0)');
%         scatter(swrPSL{a}, swrBS{a}-swrS{a}, '.')
%         xlabel('swr pct since lick');
%         ylabel('swr time since burst');
    end
end
%%

