
%{
- need to define timefilters for my two conditions:
- first condition is during lick bouts.. i have already created something
that works in the last notebook.. so i need to turn that into
timefilterfunction call
- the second condition is during immobility but outside of lick bouts..
maybe like > 500 ms away from the nearest lick.
++ ok i made getLickBout and i think/hope it works
... now how to do the analyses

lickbout VS no-lickbout(cm/s < 1)
STEP 1::: (create timefilter or intervals? for each, and verify)
===========================

INTRA AREA
HPC, MECs, MECd. Spike-Spike. "SU ExCorr":: multicellanal
HPC, MECs, MECd. Spike-LFP. "Phase locking":: singlecellanal ?
HPC, MECs, MECd. LFP-LFP. "Coherence?":: multitetanal?

INTER AREA
Ordered by Strongest demonstration (HPC-MEC)..
HPC-MECs, HPC-MECd, MECs-MECd. Spike-Spike
HPC-MECs, HPC-MECd, MECs-MECd. spike-LFP
HPC-MECs, HPC-MECX-area Spike-Spike peak ExCorr

SWR PROPAGATION

%}
clear all
create_filter = 1;
run_ff = 1;
save_ffdata = 1;

load_ffdata = 0;
stack_spikes = 0;

make_licks = 0;

pilotTimeFilter = 0;
loadEpoch = 0;
run_dfa = 0;

plotfigs = 0;
displayplots = 0;
saveplots = 1;

%% data filter params
conditions = {'lickbouts', 'nolickbouts'};
for c = 1:length(conditions)
    clear Fp F
    condition = conditions{c};
    Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
    Fp.filtfunction = 'dfa_lickBoutSpikeCorr';
    Fp.filtfunction = 'dfa_lickswrcorr';
    Fp.add_params = {'savefigs', 'wtrackdays', 'exemplar_wepochs', condition};
    %% FF
    Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
    if create_filter
        F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
            'excludetime', Fp.timefilter,'iterator',Fp.iterator,'cellpairs', Fp.cellpairfilter);
        F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    end
    if run_ff; F = runfilter(F);
        for d = 1:length(F); F(d).datafilter_params = Fp; end
    end
    if save_ffdata
        save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
            'filetail', sprintf('_%s_%s', Fp.epochEnvironment, condition));
    end
end
if load_ffdata
    lickF = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', sprintf('_%s_%s', Fp.epochEnvironment, conditions{1}));
    ctlF = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', sprintf('_%s_%s', Fp.epochEnvironment, conditions{2})); end

%% get distribution of excessdorr for each pair /area-pair /condition
% *** confirm the conditions.. i.e. spikes used for each, along with the
% lick dio plot.. cuz i think i have this backwards right now
% maybe do heatmap of all pairs lick-ctrl
% also need to show per area-pair lick-ctrl summary and tests
if plotfigs
%{
Done    % 
Figure 1:
    plot data with lick bout intervals and no bout intervals : 
    confirm conditions and get examples illustrating effect
        - plot pos, speed, licking, spiking, lfp
    - maybe just plot the datachunks with the intervals overlaid? yes
    - [plotDataChunks_20190904 : dfa_plotDataChunks]

Figure 2:     
    lick - SWR xcorr

Figure 3: lfp correlation across and within areas change during licking
    
Figure 4: spike phase locking to intra and inter area lfp
    
Figure 5: Spiking correlation across and within areas change during licking
    
    plot the xcorr +- sem for each pairing..
    
    plot box whiskers of excesscorr bout V nobout. and/or cdf/pdf distributions and tests: 
    
    plot proportion of pairs where the excess corr bout V nobout is significantly different
%}

%     
%     for ani = 1:length(lickF)
%         excorr_lick = cellfun(@(x) x.excesscorr, lickF(ani).output{1}, 'un', 1)';
%         excorr_ctrl = cellfun(@(x) x.excesscorr, ctlF(ani).output{1}, 'un', 1)';
%         [d1,~,x] = unique(lickF(ani).data{1}{1}(:,1:2), 'rows');
%         [d2,~,y] = unique(lickF(ani).data{1}{1}(:,3:4), 'rows');
%         % label each with area so that i can group by intra and inter area
%         % pairings.. should be 6 groups total..
%         
%         % how many pairs in each group are significantly correlated?
%         
%         
%         z = excorr_lick - excorr_ctrl;
%         exc_heatmap = accumarray([x(:),y(:)],z(:),[size(d1,1) size(d2,1)]);
%         imagesc(exc_heatmap)
%         colormap(jet)
%         u(x,y) = excorr_lick-excorr_ctrl;
%         heatmap(x,y,excorr_lick-excorr_ctrl)
%         
%         figure
%         subplot(2,2,1)
%         ecdf(excorr_lick(~isnan(excorr_lick)),'Bounds','on')
%         hold on
%         ecdf(excorr_ctrl,'Bounds','on')
%         %         x = linspace(min(excorr_lick),max(excorr_lick));
%         %         plot(x,evcdf(x,0,3))
%         grid on
%         title('Empirical CDF')
% 
%         hold off
%         subplot(2,2,2)
%         [f,xi] = ksdensity(excorr_lick);
%         plot(xi,f);
%         hold on
%         [f,xi] = ksdensity(excorr_ctrl);
%         plot(xi,f);
%         hold off
%         title('KDEst')
%         subplot(2,2,3)
%         
%     end
end
%%






