


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
plotETAFull = 1;
plotETAPerPC = 0;
plotTrace = 0;
plotPCdemo = 0;
pausefigs = 0;
savefigs = 1;

%% FF
Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
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
%% Concat and Create results across epochs per animal
if createSummaryData
    D = struct;
    for i = 1:length(F)
        f = [F(i).output{:}];
        D(i).animal = F(i).animal;
        D(i).etTime = F(i).output{end}.etaTime;
        D(i).swrReactETAfull = cell2mat(cellfun(@(x) x', {f.swrReactETAfull}', 'un', 0)')';
        D(i).swrReactETAfullMean = nanmean(D(i).swrReactETAfull);
%         D(i).swrReactETAfullShufs = cell2mat(cellfun(@(x) x', {f.swrReactETAfullShufs}', 'un', 0)')';
%         D(i).swrReactETAfullShufMean = nanmean(D(i).swrReactETAfullShufs);
        
        D(i).swrBurstReactETAfull = cell2mat(cellfun(@(x) x', {f.swrBurstReactETAfull}', 'un', 0)')';
        D(i).swrBurstReactETAfullMean = nanmean(D(i).swrBurstReactETAfull);
%         D(i).swrBurstReactETAfullShufs = cell2mat(cellfun(@(x) x', {f.swrBurstReactETAfullShufs}', 'un', 0)')';
%         D(i).swrBurstReactETAfullShufMean = nanmean(D(i).swrBurstReactETAfullShufs);
        
        D(i).lickReactETAfull = cell2mat(cellfun(@(x) x', {f.lickReactETAfull}', 'un', 0)')';
        D(i).lickReactETAfullMean = nanmean(D(i).lickReactETAfull);
%         D(i).lickReactETAfullShufs = cell2mat(cellfun(@(x) x', {f.lickReactETAfullShufs}', 'un', 0)')';
%         D(i).lickReactETAfullShufMean = nanmean(D(i).lickReactETAfullShufs);
        
        % per pc
        u = all([~cellfun(@isempty, {f.swrReactETAPerPC}, 'un', 1); ...
             ~cellfun(@isempty, {f.swrBurstReactETAPerPC}, 'un', 1); ...
             ~cellfun(@isempty, {f.lickReactETAPerPC}, 'un', 1)]);
         
        D(i).dayEpoch = cell2mat({f(u).idx}');
        D(i).eigValSortSig = {f(u).eigValSortSig};
        
        D(i).swrReactETAPerPC = cellfun(@(x) cell2mat(cellfun(@(y) y',x,'un',0)),{f(u).swrReactETAPerPC},'un', 0);
        D(i).swrReactETAPerPCMean = cellfun(@(x) nanmean(x,2), D(i).swrReactETAPerPC,'un', 0);
        %         D(i).swrReactETAfullShufs = cellfun(@(x) cell2mat(cellfun(@(y) y',x,'un',0)),{data.swrReactETAPerPC},'un', 0);
        %         D(i).swrReactETAfullShufMean = nanmean(D(i).swrReactETAfullShufs);
        
        D(i).swrBurstReactETAPerPC = cellfun(@(x) cell2mat(cellfun(@(y) y',x,'un',0)),{f(u).swrBurstReactETAPerPC},'un', 0);
        D(i).swrBurstReactETAPerPCMean = cellfun(@(x) nanmean(x,2), D(i).swrBurstReactETAPerPC,'un', 0);
        
        D(i).lickReactETAPerPC = cellfun(@(x) cell2mat(cellfun(@(y) y',x,'un',0)),{f(u).lickReactETAPerPC},'un', 0);
        D(i).lickReactETAPerPCMean = cellfun(@(x) nanmean(x,2), D(i).lickReactETAPerPC,'un', 0);
    end
end

%% ------------------------------plot------------------------------
%% plot perAnimal allEventGroups rxn ETA full
if plotfigs
    if plotETAFull
        figname = 'RxnFull';
        for i = 1:length(D)
            %% Fig start
            Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
                'savefigs', savefigs);
            time = D(i).etTime';
            winidx = knnsearch(time, Pp.winSE');
            time = time(winidx(1):winidx(2));
            
            %% plot all swr rxn ETA full
            Rs = D(i).swrReactETAfull(:,winidx(1):winidx(2));
            interptime = [time(1):diff(time(1:2))/Pp.interpBy:time(end)]';
            RsP = interp1(time, Rs', interptime, 'spline')';
            R = D(i).swrReactETAfullMean(:,winidx(1):winidx(2));
            RP = interp1(time, R, interptime, 'spline')';
            Rsem = nanstd(RsP)/sqrt(size(RsP,1));
%             shs = D(i).swrReactETAfullShufs(:,winidx(1):winidx(2));
%             sh = D(i).swrReactETAfullShufMean(:,winidx(1):winidx(2));
%             shssem = nanstd(shs)/sqrt(size(shs,1));            
%             plot(time, sh, 'k', 'linewidth', 1.5)
%             fill([time; flipud(time)]',[sh-shssem flipud(sh+shssem)]','k','linestyle','none', ...
%                 'facealpha', .1)
            fill([interptime; flipud(interptime)],[RP'+Rsem'; flipud(RP'-Rsem')],'k', ...
                'linestyle','none', 'facealpha', .1)
            hold on;
            f1 = plot(interptime, RP, 'k', 'linewidth', 1, 'DisplayName','SWR');
            
            %% plot swrBurst rxn ETA full
            Rs = D(i).swrBurstReactETAfull(:,winidx(1):winidx(2));
            RsP = interp1(time, Rs', interptime, 'spline')';
            R = D(i).swrBurstReactETAfullMean(:,winidx(1):winidx(2));
            RP = interp1(time, R, interptime, 'spline')';
            Rsem = nanstd(RsP)/sqrt(size(RsP,1));
            fill([interptime; flipud(interptime)],[RP'+Rsem'; flipud(RP'-Rsem')],'r', ...
                'linestyle','none', 'facealpha', .1)
            hold on;
            f2 = plot(interptime, RP, 'r', 'linewidth', 1, 'DisplayName','BurstSWR');
            
            %% plot lickBurst rxn ETA full
            Rs = D(i).lickReactETAfullMean(:,winidx(1):winidx(2));
            R = D(i).lickReactETAfullMean(:,winidx(1):winidx(2));
            RP = interp1(time, R, interptime, 'spline')';
            Rsem = nanstd(RsP)/sqrt(size(RsP,1));
            fill([interptime; flipud(interptime)],[RP'+Rsem'; flipud(RP'-Rsem')],'b', ...
                'linestyle','none', 'facealpha', .1)
            hold on;
            f3 = plot(interptime, RP, 'b', 'linewidth', 1, 'DisplayName','BurstLick');
            
            axis tight
            ylabel('rxn strength')
            xlabel('Time from lickDin (s)')
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
            legend([f1 f2 f3]);
            hold off;
            
            %% End Fig
            try
                animal = D(i).animal{3};
            catch
                animal = D(i).animal;
            end
            stit = sprintf('%s %s %s', animal, figname, Fp.env);
            setSuperAxTitle(stit);
            
            if pausefigs; pause; end
            if savefigs
                save_figure(strjoin({pconf.andef{4}, figname, animal}, '/'), stit);
            end
        end
    end
    
    %% plot perEpoch perEventGroup perPC rxn ETA
    if plotETAPerPC
        figname = 'RxnPerPC';
        for i = 1:length(D) % animal
            for iep = 1:length(D(i).swrReactETAPerPC) % epoch
                %% Fig start
                Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
                    'savefigs', savefigs);
                time = D(i).etTime';
                winidx = knnsearch(time, Pp.winSE');
                time = time(winidx(1):winidx(2));
                c = vals2rgb(1-D(i).eigValSortSig{iep}', 'parula');
                
                %% subplot epoch all swr rxn ETA perPC
                s1 = subplot(3,1,1);
                Rs = D(i).swrReactETAPerPC{iep}(winidx(1):winidx(2),:)';
                hold on;
                arrayfun(@(x) plot(time, Rs(x,:), 'color', c(x,:), 'linewidth', 1), 1:size(c,1), 'un', 0);
                axis tight
                hold off;
                line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
                ylabel('swr rxn')
                
                %% subplot epoch all swrBurst rxn ETA perPC
                s2 = subplot(3,1,2);
                Rs = D(i).swrBurstReactETAPerPC{iep}(winidx(1):winidx(2),:)';
                hold on;
                arrayfun(@(x) plot(time, Rs(x,:), 'color', c(x,:), 'linewidth', 1), 1:size(c,1), 'un', 0);
                axis tight
                hold off;
                line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
                ylabel('swrBurst rxn')
                
                %% subplot epoch lickBurst rxn ETA perPC
                s3 = subplot(3,1,3);
                Rs = D(i).lickReactETAPerPC{iep}(winidx(1):winidx(2),:)';
                hold on;
                arrayfun(@(x) plot(time, Rs(x,:), 'color', c(x,:), 'linewidth', 1), 1:size(c,1), 'un', 0);
                axis tight
                hold off;
                line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
                ylabel('lickBurst rxn')
                xlabel('Time from event (s)')
                
               %% End Fig
                allAxesInFigure = findall(gcf,'type','axes');
                linkaxes(allAxesInFigure, 'y');
                try
                    animal = D(i).animal{3};
                catch
                    animal = D(i).animal;
                end
                stit = sprintf('%s %d %d %s %s', animal, D(i).dayEpoch(iep,:), ...
                    figname, Fp.env);
                setSuperAxTitle(stit);
                
                if pausefigs; pause; end
                if savefigs
                    save_figure(strjoin({pconf.andef{4}, figname, animal}, '/'), stit);
                end
            end
        end
    end
end

%% plot spatial correspondance to the PC's perEpoch

%% plot perAnimal top 3 PC mean for swr in burst vs licks










