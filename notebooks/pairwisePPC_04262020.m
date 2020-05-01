

%{
pairwise phase/time/order consistency during lick bouts
ok.. so he approach is going to turn time into phase series and do an xcorr

how do i load like a couple cells and lick periods to just tinker before
turning things into a more formal filter function?
where is load individual epoch code blocks or eample sripts?

where is the code i've used for xcorr?

Durkewitz? book chapter 7 includes good xcorr stuff but no code..
Lauback has a xcorr code snip

how make this FF compatible, I need to make multicelldayanal
which would feed in pairs of cell data for the whole day to the dfa function

im realizing that my 'notebooks' kinda were what DFS were, except more
journal like.. so kinda like final or demo notebooks that can be archived

What is dfa_lickXCorrSpikes?
when did i use the param set lickspikexc
what about 'riptrigspiking_corr'? which has 'calcxcorrmeasures' listed as filtfunciton
what about 'dfa_suCoactiveXcorr'?? **********
what about 'mua_calcxcorrmeasures'?
what about 'dfa_perripspikingcorr'?

%}

pconf = paramconfig;
create_filter = 1;
run_ff = 0;
load_ffdata = 0;

savefigs = 0;
pausefigs = 1;
showfigs = 1;
savefigas = {'png','pdf'};

plot_phaseXcorr = 0;

%% Define Filter Params
pconf = paramconfig('Demetris'); % globals per user
Fp.animals = {'JZ4'};
eventType = 'lick';
Fp.filtfunction = 'dfa_phaseXcorr';
Fp.Label = 'wtrackLickTrigSpiking'; % used for filename/plots
Fp.params = {'exemplar_wepochs', ...
    'valid_ntrodes', ...
    'excludePriorFirstWell', ...
    'excludeAfterLastWell', ...
    'nonMU_cells', ...
    Fp.Label, ...
    Fp.filtfunction};
Fp = load_filter_params(Fp);

%% Create Filter and Run FiltFunc
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, ...
        'eegtetrodes', Fp.tetfilter, 'excludetime', Fp.timefilter, ...
        'iterator', Fp.iterator, 'cells', Fp.cellFilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F),'un',1);
    F = runfilter(F);
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.Label]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', ['_' Fp.Label]);
end

%% step 0: run createfilter to get set of data labels as var 'Fp'
% use load_filter_params to use past param sets or get ideas
% load_filter_params also sets the paths of user from animaldef into Fp.paths

%% spikes, lick bout p-events:: load epoch of data
% use load_data with 'filterframework' as source
% i need spikes and XP
% I SHOULD have been able to feed in a filter key set that createfilter SHOULD have made, like in pyphy
% instead, i have to get animals from Fp.animals, days and epochs from
% F.epochs{1} ... and ntrodes and cells from F.data{1}{1}.. but i still
% have to load all of whatever type of data before doing the filtering. 
% essentially what i want to do is what the iterators fo for a given dfa
% function.. so in essence i'm stepping into one specific dfa function call
% so maybe i instead should use an existing iterator? right now it looks
% like from the load_filter_params that singlecellanal was used for
% pairwise cell iteration.. which is stupid but is what it is for now
% so how do i take advantage of the iterator in this way to return selected data to
% the current workspace?
% the data tables that i need are listed in F.function.loadvariables

spikes = load_data('filterframework', 'spikes', Fp.animals);

% get pair spiketrains
t1 = spikes{ind(1)}{ind(2)}{ind(3)}{ind(4)}.data;
t2 = spikes{ind(1)}{ind(2)}{ind(5)}{ind(6)}.data;

% excludetime
t1inc = t1(find(~isExcluded(t1(:,1), excludetimes)),1);
t2inc = t2(find(~isExcluded(t2(:,1), excludetimes)),1);

%% transform time series into phase series


% compute xc, excorr
xc = spikexcorr(t1inc, t2inc, bin, tmax);

% compute the excess correlation at 0 lag
out.ec = excesscorr(xc.time, xc.c1vsc2, xc.nspikes1, xc.nspikes2, sw1, sw2);
out.rms = xcorrrms(xc.time, xc.c1vsc2, rmstmax, rmsmincounts);

% get the probability that both were coactive in the included intervals
n1 = nInInterval(t1inc, excludetimes, 0);
n2 = nInInterval(t2inc, excludetimes, 0);
out.probcoa = mean(logical(n1) & logical(n2));

% get the zscore for this probabliity of coactivity
out.coactivez = coactivezscore(n1, n2);
out.coactivez_numsp = coactivezscore(n1, n2, 'anyspikes', 0);

%% Notes for x-corr 

%{
% Lauback xcorr code::
edges = 0:0.001:10; % 1 ms bins
t1_binned = histc(t1, edges); % timestamps to binned
t2_binned = histc(t2, edges); % timestamps to binned
xc = xcorr(t1_binned, t2_binned); % actual crosscorrelation
%}


% mari's rip-rip xcorr dfa: dfams_ripplexcorr.m
% mari's rip-rip xcorr dfs: dfsms_ripplexcorr.m
%       look in the rip xcorr scripts for bar+error plots

% mari's spike-spike xcorr dfa: dfams_spikeCorrOutsideRip.m
%       look for spike-spike dfa for actual computation and pair scores??
%       or is this done for the day in the dfs?
% mari's spikes-spike xcorr dfs: DFScripts/dfsms_spikeCorr_outside_ripples.m

%{
mari's spike-spike xcorr methods
Spike cross-correlations between pairs of D+ or V+ NAc MSNs were
calculated in 10 ms bins at up to 0.5 s lag. Each CCH was first normalized by the square root of
the product of the number of spikes from each cell. To z-score the CCH of each cell pair, one of
the spike trains was circularly shuffled 1000 times (by a random amount up to ±half the mean
immobility period length) to create 1000 shuffled CCHs. Each real and shuffled CCH was
smoothed with a 20 ms s.d., 160 ms wide Gaussian kernel. The real cross-correlation values were
then z-scored relative to the distribution of shuffled values within each bin. We averaged the
smoothed cross-correlation z-score in a 20 ms bin around 0 (±10 ms) to get an approximate “0
%}

%% My dfa for spike-spike xcorr: 
% binsize = 10 % ms
% tmax = .5 % s
% compute xc using spikexcorr.m
% normalize by square root of the product of num of spikes from each cell
% 1000 times circ shuf spike trains up to half of lick-bout period
% smooth real and shuffled with 20 ms s.d., 160 ms wide gaussian kernel
% zscore real relative to distribution of shuffle vals within each time bin
% DON'T:: average the result in a 20 ms bin around 0 to get approximate "0 lag".. 
% INSTEAD:: for me, get the max xcorr within half the lower length of p-event
% period.. so like 25 or 30 ms 

%% plot xcorr of each cell pair
if plot_phaseXcorr
    figname = 'phaseXcorr';
    Pp = load_plotting_params({'defaults', figname}); % load params
    for a = 1:length(F) % per animal
%         animal = F(a).animal{3};
%         pairs = F.pairs;
        for i = 1:size(pairs,1) % per pairwise phase xcorr
            ifig = init_plot(showfigs, Pp.position); % init fig
            sf = subaxis(1,1,1,Pp.posparams{:});
            sf.Tag = 'xcorr';
            % plot bar and errorbar binned time smoothed zscored xcorr 
            % Z-score 0 reflects mean of shuffles. Error bars indicate s.e.m.
            
            % superfig
            stit = sprintf('%s %s %s %s', figname, animal, Fp.env, ....
                strjoin(Fp.areas{:,i}));
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure([pconf.andef{4} '/' figname '/' animal], ...
                    stit, 'savefigas');
            end
        end
    end
end
%% plot cumulative peak xcorr for pairs vs shuff
% mari's uses one-tailed Wilcoxon rank-sum test
%% plot fraction of pairs sig xcorr
% mari's uses a z-test for proportions
%% plot xcorr matrix of all pairs