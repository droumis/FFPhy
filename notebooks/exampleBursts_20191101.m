


%{
plot example of a swr

plot example of lick burst with swr's

also today:
- swr (burst, nonburst) reactivation
- su swr mod (burst, nonburst), lick su mod
- wet lick burst vs dry lick burst

i need to confirm that the swr, licks, im using for all the analysis are
correctly filtered for.. via looking at the example traces..

i need to confirm the reactivation trace result along with the example
traces
%}

pconf = paramconfig;
create_filter = 1;
run_ff = 1;
load_ffdata = 0;
createSummaryData = 1;

plotfigs = 1;
plotETAFull = 0;
plotETAPerPC = 1;
plotTrace = 0;
plotPCdemo = 0;
pausefigs = 0;
savefigs = 1;

%% FF
Fp.animals = {'JZ1'};
Fp.filtfunction = 'dfa_reactivationPSTH';
Fp.params = {'ripples', 'ca1SU', 'wtrackdays', 'excludePriorFirstWell', Fp.filtfunction};
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
%% create summary data
if createSummaryData
    D = struct;
    for i = 1:length(F)
        f = [F(i).output{:}];
        D(i).animal = F(i).animal;
        D(i).etTime = F(i).output{end}.etaTime;
        D(i).swrReactETAfull = cell2mat(cellfun(@(x) x', {f.swrReactETAfull}', 'un', 0)')';
        D(i).swrReactETAfullMean = nanmean(D(i).swrReactETAfull);
        D(i).swrReactETAfullShufs = cell2mat(cellfun(@(x) x', {f.swrReactETAfullShufs}', 'un', 0)')';
        D(i).swrReactETAfullShufMean = nanmean(D(i).swrReactETAfullShufs);
        
        D(i).swrBurstReactETAfull = cell2mat(cellfun(@(x) x', {f.swrBurstReactETAfull}', 'un', 0)')';
        D(i).swrBurstReactETAfullMean = nanmean(D(i).swrBurstReactETAfull);
        D(i).swrBurstReactETAfullShufs = cell2mat(cellfun(@(x) x', {f.swrBurstReactETAfullShufs}', 'un', 0)')';
        D(i).swrBurstReactETAfullShufMean = nanmean(D(i).swrBurstReactETAfullShufs);
        
        D(i).lickReactETAfull = cell2mat(cellfun(@(x) x', {f.lickReactETAfull}', 'un', 0)')';
        D(i).lickReactETAfullMean = nanmean(D(i).lickReactETAfull);
        D(i).lickReactETAfullShufs = cell2mat(cellfun(@(x) x', {f.lickReactETAfullShufs}', 'un', 0)')';
        D(i).lickReactETAfullShufMean = nanmean(D(i).lickReactETAfullShufs);
        
        % per pc
        D(i).eigValSortSig = {f.eigValSortSig};
        
        D(i).swrReactETAPerPC = cellfun(@(x) cell2mat(cellfun(@(y) y',x,'un',0)),{f.swrReactETAPerPC},'un', 0);
        D(i).swrReactETAPerPCMean = cellfun(@(x) nanmean(x,2), D(i).swrReactETAPerPC,'un', 0);
        %         D(i).swrReactETAfullShufs = cellfun(@(x) cell2mat(cellfun(@(y) y',x,'un',0)),{data.swrReactETAPerPC},'un', 0);
        %         D(i).swrReactETAfullShufMean = nanmean(D(i).swrReactETAfullShufs);
        D(i).swrBurstReactETAPerPC = cellfun(@(x) cell2mat(cellfun(@(y) y',x,'un',0)),{f.swrBurstReactETAPerPC},'un', 0);
        D(i).swrBurstReactETAPerPCMean = cellfun(@(x) nanmean(x,2), D(i).swrBurstReactETAPerPC,'un', 0);
        
        D(i).lickReactETAPerPC = cellfun(@(x) cell2mat(cellfun(@(y) y',x,'un',0)),{f.lickReactETAPerPC},'un', 0);
        D(i).lickReactETAPerPCMean = cellfun(@(x) nanmean(x,2), D(i).lickReactETAPerPC,'un', 0);
    end
end



%% plot
if plotfigs
    if plotETAFull
        %% plot Full ETA w sem errorbars per animal
        figname = 'Reactivation Full';
        for i = 1:length(D)
            Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
                'savefigs', savefigs);
            
            time = D(i).etTime';
            winidx = knnsearch(time, Pp.winSE');
            time = time(winidx(1):winidx(2));
            
            subplot(1,1,1)
            Rs = D(i).swrReactETAfull(:,winidx(1):winidx(2));
            R = D(i).swrReactETAfullMean(:,winidx(1):winidx(2));
            shs = D(i).swrReactETAfullShufs(:,winidx(1):winidx(2));
            sh = D(i).swrReactETAfullShufMean(:,winidx(1):winidx(2));
            Rsem = nanstd(Rs)/sqrt(size(Rs,1));
            shssem = nanstd(shs)/sqrt(size(shs,1));
            
            plot(time, R, 'g', 'linewidth', 1.5)
            hold on;
            plot(time, sh, 'k', 'linewidth', 1.5)
            fill([time; flipud(time)],[R'+Rsem'; flipud(R'-Rsem')],'g','linestyle','none', ...
                'facealpha', .1)
            fill([time; flipud(time)]',[sh-shssem flipud(sh+shssem)]','k','linestyle','none', ...
                'facealpha', .1)
            ylabel('rxn strength')
            
            % plot swr burst ETA with ebars, shuff
            Rs = D(i).swrBurstReactETAfull(:,winidx(1):winidx(2));
            R = D(i).swrBurstReactETAfullMean(:,winidx(1):winidx(2));
            shs = D(i).swrBurstReactETAfullShufs(:,winidx(1):winidx(2));
            sh = D(i).swrBurstReactETAfullShufMean(:,winidx(1):winidx(2));
            Rsem = nanstd(Rs)/sqrt(size(Rs,1));
            shssem = nanstd(shs)/sqrt(size(shs,1));
            
            plot(time, R, 'b', 'linewidth', 1.5)
            hold on;
            plot(time, sh, 'k', 'linewidth', 1.5)
            fill([time; flipud(time)],[R'+Rsem'; flipud(R'-Rsem')],'b','linestyle','none', ...
                'facealpha', .1)
            fill([time; flipud(time)]',[sh-shssem flipud(sh+shssem)]','k','linestyle','none', ...
                'facealpha', .1)
            
            % plot lick ETA with ebars, shuff
            Rs = D(i).lickReactETAfull(:,winidx(1):winidx(2));
            R = D(i).lickReactETAfullMean(:,winidx(1):winidx(2));
            shs = D(i).lickReactETAfullShufs(:,winidx(1):winidx(2));
            sh = D(i).lickReactETAfullShufMean(:,winidx(1):winidx(2));
            Rsem = nanstd(Rs)/sqrt(size(Rs,1));
            shssem = nanstd(shs)/sqrt(size(shs,1));
            
            plot(time, R, 'm', 'linewidth', 1.5)
            hold on;
            plot(time, sh, 'k', 'linewidth', 1.5)
            fill([time; flipud(time)],[R'+Rsem'; flipud(R'-Rsem')],'m','linestyle','none', ...
                'facealpha', .1)
            fill([time; flipud(time)]',[sh-shssem flipud(sh+shssem)]','k','linestyle','none', ...
                'facealpha', .1)
            axis tight
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
            ylabel('rxn strength')
            xlabel('Time from lickDin (s)')
            hold off;
            
            %%
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
    %% plot per PC ETA
    if plotETAPerPC
        figname = 'Reactivation PerPC';
        for i = 1:length(D) % animal
            for iep = 1:length(D(i).swrReactETAPerPC) % epoch
                Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
                    'savefigs', savefigs);
                time = D(i).etTime';
                winidx = knnsearch(time, Pp.winSE');
                time = time(winidx(1):winidx(2));
                c = vals2rgb(1-D(i).eigValSortSig{iep}', 'parula');
                %%
                s1 = subplot(3,1,1);
                Rs = D(i).swrReactETAPerPC{iep}(winidx(1):winidx(2),:)';
                hold on;
                arrayfun(@(x) plot(time, Rs(x,:), 'color', c(x,:), 'linewidth', 1), 1:size(c,1), 'un', 0);
                axis tight
                ylabel('swr rxn')
                
                %% plot swr burst ETA with ebars, shuff
                s2 = subplot(3,1,2);
                Rs = D(i).swrBurstReactETAPerPC{iep}(winidx(1):winidx(2),:)';
                hold on;
                arrayfun(@(x) plot(time, Rs(x,:), 'color', c(x,:), 'linewidth', 1), 1:size(c,1), 'un', 0);
                axis tight
                ylabel('swrBurst rxn')
                
                %% plot lick ETA with ebars, shuff
                s3 = subplot(3,1,3);
                Rs = D(i).lickReactETAPerPC{iep}(winidx(1):winidx(2),:)';
                hold on;
                arrayfun(@(x) plot(time, Rs(x,:), 'color', c(x,:), 'linewidth', 1), 1:size(c,1), 'un', 0);
                axis tight
                ylabel('lickBurst rxn')
                xlabel('Time from event (s)')
                line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
                hold off;
                %%
                allAxesInFigure = findall(gcf,'type','axes');
                linkaxes(allAxesInFigure, 'y');
                try
                    stit = [figname ' ' D(i).animal{3} ' ' Fp.env];
                catch
                    stit = [figname ' ' D(i).animal ' ' Fp.env];
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
    end
end








