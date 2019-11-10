%{

So far: 
    - dfa_eventTrigSpiking (collected psth and eta; singlecellanal) 
        - combineEpochs
        - calcSUmod

    - dfa_lickXCorrSpikes (xcorr, excorr, phasemod; singleDayCellAnal) 

need to combine the sumod/psth code with the lickXcorrSpikes code.. 
- or do i?
- ideally dfa_eventTrigSpiking would call calcSUmod within the FF, but for
that i would need to combine the epochs via singleDayCellAnal and therefore
refactor dfa_eventTrigSpiking and calcSUmod a bit.. 

idk what to do.. i guess i just want some examples of lick triggered
spiking from SU, and then a population summary of the modulation score for
the population of SU.. and also a ETA of all the SU cells.. which i guess
could either be the ETA from dfa_eventTrigSpiking or the xcorr from
dfa_lickXCorrSpikes... 
the 'mod' score could be 
    - par: Raleigh, nonpar: excorr, phasemod, sumod..

1. remake the lick triggered rasters with dfa_eventTrigSpiking via singleDayCellAnal
2. 
%}

create_filter = 0;
run_ff = 1;
load_ffdata = 0;

pausefigs = 0;
savefigs = 1;
Fp.Label = 'lickTrigSUmod';
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_eventTrigSpiking';
Fp.params = {Fp.Label, 'savefigs', 'wtrackdays', 'valid_ntrodes', Fp.filtfunction, 'nonMU_cells'};
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
        'filetail', sprintf('_%s_%s', Fp.env))
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', sprintf('_%s', Fp.epochEnvironment));
end




