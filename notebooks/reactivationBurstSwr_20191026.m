

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
pconf = paramconfig;
create_filter = 0;
run_ff = 0;
load_ffdata = 0;
createSummaryData = 1;

plotfigs = 1;
plotETA = 1;
plotTrace = 0;
plotPCdemo = 0;
pausefigs = 0;
savefigs = 1;
%% FF
Fp.animals = {'D10'};
Fp.filtfunction = 'dfa_reactivationPSTH';
Fp.params = {'>4cm/s', 'ca1SU', 'wtrackdays', 'excludePriorFirstWell', Fp.filtfunction};
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
    D = struct;
    for i = 1:length(F)
        data = [F(i).output{:}];
        D(i).animal = F(i).animal;
        D(i).etTime = F(i).output{end}.reactPSTHtime;
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
    if plotETA
        %% plot Full ETA w sem errorbars
        figname = 'ReactivationStrength';
        for i = 1:length(D)
            Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
                'savefigs', savefigs);
%             plot swr ETA with ebars, shuff
            time = D(i).etTime';
            subplot(1,1,1)
            winidx = knnsearch(time, Pp.winSE');
            time = time(winidx(1):winidx(2));
            Rs = D(i).swrReactETAfull(:,winidx(1):winidx(2));
            R = D(i).swrReactETAfullMean(:,winidx(1):winidx(2));
            shs = D(i).swrReactETAfullShufs(:,winidx(1):winidx(2));
            sh = D(i).swrReactETAfullShufMean(:,winidx(1):winidx(2));
            Rsem = nanstd(Rs)/sqrt(size(Rs,1));
            shssem = nanstd(shs)/sqrt(size(shs,1));
            plot(time, R, 'g', 'linewidth', 1.5)
%             errorbar(time, R, Rsem, 'b','linewidth',0.5,'CapSize',0)
            hold on;
            plot(time, sh, 'k', 'linewidth', 1.5)
%             errorbar(time, sh, shssem,'k', 'linewidth',0.5, 'CapSize',0)
            fill([time; flipud(time)],[R'+Rsem'; flipud(R'-Rsem')],'g','linestyle','none', ...
                'facealpha', .1)
            fill([time; flipud(time)]',[sh-shssem flipud(sh+shssem)]','k','linestyle','none', ...
                'facealpha', .1)
%             axis tight
%             line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
            ylabel('rxn strength')
%             xticks([]);
%             title('SWR')
%             hold off;
            
            %% plot swr burst ETA with ebars, shuff
            subplot(1,1,1)
            Rs = D(i).swrBurstReactETAfull(:,winidx(1):winidx(2));
            R = D(i).swrBurstReactETAfullMean(:,winidx(1):winidx(2));
            shs = D(i).swrBurstReactETAfullShufs(:,winidx(1):winidx(2));
            sh = D(i).swrBurstReactETAfullShufMean(:,winidx(1):winidx(2));
            Rsem = nanstd(Rs)/sqrt(size(Rs,1));
            shssem = nanstd(shs)/sqrt(size(shs,1));
            plot(time, R, 'b', 'linewidth', 1.5)
%             errorbar(time, R, Rsem, 'b','linewidth',0.5,'CapSize',0)
            hold on;
            plot(time, sh, 'k', 'linewidth', 1.5)
%             errorbar(time, sh, shssem,'k', 'linewidth',0.5, 'CapSize',0)
            fill([time; flipud(time)],[R'+Rsem'; flipud(R'-Rsem')],'b','linestyle','none', ...
                'facealpha', .1)
            fill([time; flipud(time)]',[sh-shssem flipud(sh+shssem)]','k','linestyle','none', ...
                'facealpha', .1)
            ylabel('rxn strength')
%             xticks([]);
            axis tight
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
            ylabel('rxn strength')
            xlabel('Time from SWR (s)')
            title('rxn SWR ETA')
%             hold off;
            %% plot lick ETA with ebars, shuff
            subplot(1,1,1)
            Rs = D(i).lickReactETAfull(:,winidx(1):winidx(2));
            R = D(i).lickReactETAfullMean(:,winidx(1):winidx(2));
            shs = D(i).lickReactETAfullShufs(:,winidx(1):winidx(2));
            sh = D(i).lickReactETAfullShufMean(:,winidx(1):winidx(2));
            Rsem = nanstd(Rs)/sqrt(size(Rs,1));
            shssem = nanstd(shs)/sqrt(size(shs,1));
            plot(time, R, 'm', 'linewidth', 1.5)
%             errorbar(time, R, Rsem, 'b','linewidth',0.5,'CapSize',0)
            hold on;
            plot(time, sh, 'k', 'linewidth', 1.5)
%             errorbar(time, sh, shssem,'k', 'linewidth',0.5, 'CapSize',0)
            fill([time; flipud(time)],[R'+Rsem'; flipud(R'-Rsem')],'m','linestyle','none', ...
                'facealpha', .1)
            fill([time; flipud(time)]',[sh-shssem flipud(sh+shssem)]','k','linestyle','none', ...
                'facealpha', .1)
            axis tight
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
            ylabel('rxn strength')
            xlabel('Time from lickDin (s)')
            title('rxn Lick ETA')
%             legend({'allSWR', 'allSWRshuff', 'burstSWR', 'burstSWRshuff', 'lick', 'lickshuf'});
            hold off;
            %%
            allAxesInFigure = findall(gcf,'type','axes');
            linkaxes(allAxesInFigure, 'y');
            try
                stit = [figname D(i).animal{3} Fp.env];
            catch
                stit = [figname D(i).animal Fp.env];
            end
            setSuperAxTitle(stit);
            if pausefigs
                pause;
            end
            if savefigs
                save_figure([pconf.andef{4} '/' figname], stit);
            end
        end
    end
    if plotTrace
        %% plot time snippet of spikes + zbinSpikes traces + reactivation result (full,pc)
        exWin = [];
        % plot spikes
        
        % plot full reactivation trace
        
        % plot lick times
        
        % plot swr times
    end
    %% plot per PC method Demo
    % plot the steps and examples of the main parts
    % templateCC Eigvec, val, sig PC's, match z spikes, per pc react
    if plotPCdemo
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











