

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
create_filter = 1;
run_ff = 1;
load_ffdata = 0;
createSummaryData = 0;

plotfigs = 0;
pausefigs = 1;
savefigs = 0;
%% FF
Fp.animals = {'D10', 'D12', 'JZ1', 'JZ2', 'JZ4'}; %, 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_reactivationPSTH';
Fp.params = {'>4cm/s', 'ca1SU', 'wtrackdays', 'excludePriorFirstWell', ...
    Fp.filtfunction}; %'exemplar_wepochs', 
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
%% copy ntrode info to cell info for the rest of the animals



%% create summary data
if createSummaryData
    for i = 1:length(F)
        data = [F(i).output{:}];
        D(i).animal = F(i).animal;
        D(i).etTime = F(i).output{1}.reactPSTHtime;
        D(i).swrReactETAfull = cell2mat(cellfun(@(x) x', {data.swrReactETAfull}', 'un', 0)')';
        D(i).swrReactETAfullMean = nanmean(D(i).swrReactETAfull);
        D(i).swrReactETAfullShufs = cell2mat(cellfun(@(x) x', {data.swrReactETAfullShufs}', 'un', 0)')';
        D(i).swrReactETAfullShufMean = nanmean(D(i).swrReactETAfullShufs);
        
        D(i).swrBurstReactETAfull = cell2mat(cellfun(@(x) x', {data.swrBurstReactETAfull}', 'un', 0)')';
        D(i).swrBurstReactETAfullMean = nanmean(D(i).swrBurstReactETAfull);
        D(i).swrBurstReactETAfullShufs = cell2mat(cellfun(@(x) x', {data.swrBurstReactETAfullShufs}', 'un', 0)')';
        D(i).swrBurstReactETAfullShufMean = nanmean(D(i).swrBurstReactETAfullShufs);
        
        D(i).lickReactETAfull = cell2mat(cellfun(@(x) x', {data.lickReactETAfull}', 'un', 0)')';
        D(i).lickReactETAfullMean = nanmean(D(i).lickReactETAfull);
        D(i).lickReactETAfullShufs = cell2mat(cellfun(@(x) x', {data.lickReactETAfullShufs}', 'un', 0)')';
        D(i).lickReactETAfullShufMean = nanmean(D(i).lickReactETAfullShufs);
    end
end

%% plot
if plotfigs
    %% plot Full ETA w sem errorbars
    figname = 'ReactivationStrength';
    for i = 1:length(D)
        Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
            'savefigs', savefigs);
        % plot swr ETA with ebars, shuff
        time = D(i).etTime;
        subplot(3,1,1)
        plot(time, D(i).swrReactETAfullMean,'k','linewidth',3);
        hold on;
        plot(time, D(i).swrReactETAfullShufMean,'r','linewidth',3);
        errorbar(time, D(i).swrReactETAfullMean, nanstd(D(i).swrReactETAfull)/...
            sqrt(size(D(i).swrReactETAfull,1)), 'k','linewidth',0.5)
        errorbar(time, D(i).swrReactETAfullShufMean, nanstd(D(i).swrReactETAfullShufs)/...
            sqrt(size(D(i).swrReactETAfullShufs,1)), 'r','linewidth',0.5)
        ylabel('rxn strength')
        xticks([]);
%         xlabel('Time (s)')
        title('SWR')
        hold off;

        % plot swr burst ETA with ebars, shuff
        subplot(3,1,2)
        plot(time, D(i).swrBurstReactETAfullMean,'k','linewidth',3);
        hold on;
        plot(time, D(i).swrBurstReactETAfullShufMean,'r','linewidth',3);
        errorbar(time, D(i).swrBurstReactETAfullMean,nanstd(D(i).swrBurstReactETAfull)/...
            sqrt(size(D(i). swrBurstReactETAfull,1)), 'k','linewidth',0.5)
        errorbar(time, D(i).swrBurstReactETAfullShufMean, ...
            nanstd(D(i).swrBurstReactETAfullShufs)/...
            sqrt(size(D(i).swrBurstReactETAfullShufs,1)),'r','linewidth',0.5)
        ylabel('rxn strength')
%         xlabel('Time (s)')
        xticks([]);
        title('SWR in Burst')
        hold off;

        % plot lick ETA with ebars, shuff
        subplot(3,1,3)
        plot(time, D(i).lickReactETAfullMean,'k','linewidth',3);
        hold on;
        plot(time,D(i).lickReactETAfullShufMean,'r','linewidth',3);
        errorbar(time,D(i).lickReactETAfullMean,nanstd(D(i).lickReactETAfull)/...
            sqrt(size(D(i).lickReactETAfull,1)),'k','linewidth',0.5)
        errorbar(time,D(i).lickReactETAfullShufMean,nanstd(D(i).lickReactETAfullShufs)/...
            sqrt(size(D(i).lickReactETAfullShufs,1)),'r','linewidth',0.5)
        ylabel('rxn strength')
        xlabel('Time (s)')
        title('lick in Burst')
        hold off;
        %%
        stit = [figname D(i).animal{3} Fp.env];
        setSuperAxTitle(stit);
        if pausefigs
            pause;
        end
        if savefigs
            save_figure([pconf.andef{4} '/' figname '/'], stit);
        end
    end
    
    %% plot time snippet of spikes + zbinSpikes traces + reactivation result (full,pc)
    exWin = [];
    % plot spikes in time window
    
    % plot zbinspikes traces in time window
    % plot like lfp
    
    % plot reactivation result in time window
    % reactPerPC (1 trace per sig pc)
    % reactFull
    
    %% plot per PC method Demo
    % plot the steps and examples of the main parts
    % templateCC Eigvec, val, sig PC's, match z spikes, per pc react
    if 0
        figure
        subplot(2,3,1); imagesc(zCCtemplateNoDiag); ax = gca; ax.YDir = 'normal'; title('zCCtemplateNoDiag')
        ylabel('cell')
        subplot(2,3,2); imagesc(eigVal); ax = gca; ax.YDir = 'normal'; title('eigval')
        subplot(2,3,3); imagesc(eigVec); ax = gca; ax.YDir = 'normal'; title('eigvec')
        for iEig = 1:numEigValSig
            subplot(2,3,iEig+3); imagesc(squeeze(PCsNoDiag(iEig,:,:))); ax = gca; ax.YDir = 'normal';
            title(sprintf('zCCtemplate SigPC%dnoDiag', iEig))
        end
    end
    
    %% Plot reactStrength swr-psth, lick-psth, burst-swr-psth
    % plot burst swr-psth, burst lick-psth, nonburst swr-psth
    % underneath heatrasters, plot ETA, stdfill, shuffle
    
    %%
    % plot the spatial representation of the significant PC's, a la Peyrache
end

















