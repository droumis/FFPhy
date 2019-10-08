%{
 SWR-trig SU modulation for 2 conditions:
   - SWRs within lick bout intervals
   - SWRs at well but more than X s away from lick vout intervals

singlecell anal
- 2 timefilters
- how did i do the lick interval thing already.. i think i've saved lick
events so now i can pass them in as a datatype.. in order to get the lick
bout intervals i think i call getLickBout.m for a given epoch..
- i prob don't want to have to do that for all single cells since all
within each epoch will use the same lick bout intervals.. so either i could
save the intervals as another datatype and use singlecellanal..
- or use single epoch anal and nest the single cell loop myself within the
dffunction
- /home/droumis/Src/Matlab/filterframework_dr/notebooks/boutSpikeLFPLocking_20190921.m
- that was how i attempted to do the two timefilter conditions

- fuck... since i want to combine across epochs.. i need to just use the ff
to collect the swrtrig SU spiking, and then in the script i need to compute
mod, plot figs, etc

%}
pconf = paramconfig;
create_filter = 0;
run_ff = 0;
save_ffdata = 0;
load_ffdata =0;
combineEpochs =0;
saveCombinedEpochs = 0;
loadCombinedEpochs = 1;
plotfigs = 0;
% conditions = {'lickbouts', 'nolickbouts'};
% for c = 1:length(conditions)
%     clear Fp F
%     condition = conditions{c};
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
%     Fp.filtfunction = 'dfa_lickBoutSpikeCorr';
%     Fp.filtfunction = 'dfa_lickswrcorr';
Fp.filtfunction = 'dfa_riptrigspiking';
Fp.params = {'savefigs', 'wtrackdays', 'exemplar_wepochs','valid_ntrodes',...
    Fp.filtfunction, 'nonMU_cells'};
Fp = load_filter_params(Fp);
%% FF
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', Fp.tetfilter, ...
        'excludetime', Fp.timefilter,'iterator',Fp.iterator, 'cells', Fp.cellfilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff; F = runfilter(F);
    for d = 1:length(F); F(d).datafilter_params = Fp;
    end
end
if save_ffdata
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', sprintf('_%s_%s', Fp.epochEnvironment, condition))
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', sprintf('_%s', Fp.epochEnvironment));
end
% ---------------- combine epochs ---------------------------------------------
if combineEpochs
    if exist('F', 'var')
        ppF = combine_epochs(F, Fp, saveCombinedEpochs, Fp.paths);
    else
        error('create or load data filter output to combine epochs \n')
    end
end
% ---------------- Load combined epochs ---------------------------------------------------
if loadCombinedEpochs
    ppF = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals,...
        'filetail', '_combEps');
end
%% stack spikes
% stack_riptrigspiking isn't what i want to use.. it doesn't combine
% clusters across epochs.. 
% if stack_spikes
%     spikestack = stack_riptrigspiking(F);
% end
% %%
% if plotfigs
%     ian = 1;
%     day = 6;
%     epoch = 2;
%     ntrode = 2;
%     try
%         suspikes = squeeze(spikestack(ian).suspikes{day}{epoch}(:,:,ntrode))';
%     catch
%         suspikes = [];
%     end
%     [xx, yy] = find(suspikes);
%     s1 = scatterhist(xx/1000-1.001,yy, 'Kernel', 'on', 'bandwidth', ...
%         [.02; 1], 'location','SouthEast', 'Direction', 'in', 'Marker', ...
%         '+', 'color', 'k', 'MarkerSize', 2);
%     % s1(3).Position(3) = .08;
%     % s1(2).Position(4) = .08;
%     axis tight
%     ylabel('licknum','FontSize',8,'FontWeight','bold', 'FontName','Arial')
%     xlim([-Fp.win(1) Fp.win(2)]); xticks([-Fp.win:.2:Fp.win]);
%     xlabel('time s','FontSize',8,'FontWeight','bold', 'FontName','Arial');
%     line([0 0],ylim, 'color','red', 'linewidth', 1)
% end




