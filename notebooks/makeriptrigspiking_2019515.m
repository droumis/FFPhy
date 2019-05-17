
% also looking at dfs_riptrigspiking_20180423 and notebook_templaye

% make theriptrig spiking for D12, JZ2, JZ4
% i also want to make the rip trig spiking for sleep and openfield epochs
% aftewards

Fp = load_filter_params({'wtrack','ripples','nonref_ntrodes',...
    '>100spikes_cells','riptrigspiking'});
Fp.animals = {'D12'}; % don't use D!2 day 7 until the sorting is done
Fp.days = [1:6];

runFilterFramework = 1;
saveFilterOutput = 1;
loadFilterOutput = 0;
plotstuff = 0;

%% ---------------- Paths ---------------------------------------------------
paths = make_paths(Fp.filtfunction, Fp.epochEnvironment);
%% ---------------- Run FIlter -----------------------------------------------
if runFilterFramework  
    F = createfilter('animal',Fp.animals,'days', Fp.days,'epochs', ... 
        Fp.epochfilter, 'cells',Fp.cellfilter, 'excludetime', Fp.timefilter, ...
        'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, {'spikes', ...
        Fp.eventDataLabel ,'pos','task'}, 'TF',Fp.TF,'window',Fp.window, ...
        'binsize',Fp.binsize,'frbinsize',Fp.frbinsize,'minthresh', ...
        Fp.minstdthresh,'maxvelocity',Fp.maxvelocity,'minvelocity', ...
        Fp.minvelocity, 'consensus_numtets',Fp.consensus_numtets,'welldist', ...
        Fp.welldist);
    F = runfilter(F);
    for a = 1:length(F) % save filter detailes along with results
        F(a).datafilter_params = Fp;
    end
end
%% ---------------- Save Filter Output ----------------------------------------
if saveFilterOutput
    save_filter_output(F, paths.filtOutputDirectory, paths.filenamesave)
end
%% ---------------- Load Filter Output ----------------------------------------
if loadFilterOutput
    F = load_filter_output(paths.filtOutputDirectory, paths.filenamesave);
end