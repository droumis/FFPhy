


pconf = paramconfig;
create_filter = 1;
run_ff = 1;
load_ffdata = 0;
% combineEpochs = 0;

plot_XPmodSWR_per_epoch = 0;
plot_XPmodSWR_per_animal = 0;

plotfigs = 0;
pausefigs = 1;
savefigs = 1;
savefigas = {'png', 'eps'};

%% data filter params
Fp.animals = {'D10'}; %, 'D13', 'JZ1', 'JZ4'};
Fp.filtfunction = 'dfa_lickswrcorr';
Fp.params = {'wtrackdays','excludePriorFirstWell','excludeAfterLastWell', Fp.filtfunction};

%% FF
Fp = load_filter_params(Fp);
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'iterator',Fp.iterator);
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
        'filetail', sprintf('_%s', Fp.env));
end

%% plot per epoch
if plotfigs
    if plot_XPmodSWR_per_epoch == 1
        for ani = 1:length(F)
            for e = 1:length(F(ani).output)
                [~, fname,~] = fileparts(mfilename('fullpath'));
                outdir = sprintf('%s/%s/', pconf.andef{4},fname,animal);
                Pp=load_plotting_params({'defaults','dfa_lickswrcorr'}, 'pausefigs', pausefigs, ...
                    'savefigs', savefigs);
                
                %% Xcorr norm smooth, shuff mean/std
                sf1 = subaxis(2,2,1,Pp.posparams{:});
                sf1.Tag = 'xcorr';
                
                % shuffled xcorr with std ghost trail
                xmsh = mean(cell2mat(out.smtendhxcShuf'));
                xstdsh = std(cell2mat(out.smthxcShuf')); %/size(cell2mat(out.smthxcShuf'),1);
                plot(out.xc.time, xmsh, 'color', [0 0 1 .2], 'linewidth', 1);
                hold on;
                fill([out.xc.time'; flipud(out.time')],[xmsh'-xstdsh';flipud(xmsh'+xstdsh')],'b', 'linestyle', ...
                    'none', 'facealpha', .2);
                % xcorr norm
                bar(out.xc.time, out.normxc, 'k', 'facealpha', .2, 'edgealpha', 0)
                % xcorr norm smooth
                plot(out.time, out.smthxc, 'k')
                line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
                
                ylabel('xcorr');
                xlabel('time from lick s');
                %     title('xcorr
                hold off;
                %% excorr over shuff excorr cdf distr
                % relative swr from last lick
                sf2 = subaxis(2,2,2);
                histogram(cell2mat(out.excorrShuf), 60,'Normalization','probability','edgealpha', 0, 'facecolor', 'k');
                excsort = sort(cell2mat(out.excorrShuf));
                idxsig = round(sigpct*length(out.excorrShuf));
                line([excsort(idxsig) excsort(idxsig)], ylim, 'color', [0 0 0 .8], 'linestyle', '--');
                hold on;
                line([out.excorr out.excorr], ylim, 'color', 'r');
                excp = 1-sum(out.excorr>cell2mat(out.excorrShuf))/length(out.excorrShuf);
                title(sprintf('excorr %.03f p%.03f', out.excorr, excp));
                ylabel('probability')
                xlabel('excess corr')
                axis tight
                hold off;
                %% polar distr phase clustering, swrLickPhase, meanMRVmag
                sf3 = subaxis(2,2,3);
                %     a = polarhistogram(swrLickPhase, 16, 'Normalization', 'pdf', 'edgealpha', 0,...
                %         'facealpha', .5);
                polarplot([zeros(size(out.swrLickPhase,1),1) out.swrLickPhase]', ...
                    repmat([0 1],size(out.swrLickPhase,1),1)', 'color', [0 0 0 .4], 'linewidth', 4);
                hold on
                polarplot([0; out.vecang], [0; out.meanMRVmag], 'color', [1 0 .3], 'linewidth', 4)
                grid off
                rticks([])
                thetaticks([])
                title('swr ILI-phase')
                hold off
                axis tight
                %% phase mod
                sf4 = subaxis(2,2,4);
                histogram(cell2mat(out.phasemodShuf), 100, 'Normalization', 'probability', 'edgealpha', 0, 'facecolor', 'k');
                hold on;
                mrvsort = sort(cell2mat(out.phasemodShuf));
                idxsig = round(sigpct*length(out.phasemodShuf));
                line([mrvsort(idxsig) mrvsort(idxsig)], ylim, 'color', [0 0 0 .8], 'linestyle', '--');
                hold on
                line([out.phasemod out.phasemod], ylim, 'color', 'r');
                modp = 1-sum(out.phasemod>cell2mat(out.phasemodShuf))/length(out.phasemodShuf);
                title(sprintf('logMRVmag %.03f p%.03f Rpval%.03f', out.phasemod, modp, pval));
                ylabel('probability')
                xlabel('log(Rayleigh Z)')
                axis tight
                hold off
                
                %% ---- super axis -----
                sprtit = sprintf('%s %d %d %s', animal, day, epoch, fname(5:end));
                setSuperAxTitle(sprtit)
                % ---- pause, save figs ----
                if pausefigs; pause; end
                if savefigs; save_figure(outdir, sprtit, 'savefigas', savefigas); end
            end
        end
    end
end

%% plot per animal
if plot_XPmodSWR_per_animal
    for ani = 1:length(F)
        animal = F(ani).animal{3};
        if plotfigs
            Pp=load_plotting_params({'defaults','dfa_lickswrcorr'}, 'pausefigs', pausefigs, ...
                'savefigs', savefigs);
            
            %% Xcorr norm smooth, shuff mean/std
            sf1 = subaxis(2,2,1,Pp.posparams{:});
            sf1.Tag = 'xcorr';
            % shuffled xcorr mean.std ghost trail
            a = [F(ani).output{:}];
            nonempty = find(~cell2mat(cellfun(@isempty, {a.time}, 'un', 0)));
            time = cell2mat({a(nonempty(1)).time}');
            % xcorr norm smooth
            smxcmean = nanmean(cell2mat({a.smthxc}'),1);
            smxcstd = nanstd(cell2mat({a.smthxc}'),[],1);
            plot(time, smxcmean, 'k')
            hold on;
            fill([time'; flipud(time')],[smxcmean'-smxcstd';flipud(smxcmean'+smxcstd')],'k',...
                'linestyle','none', 'facealpha', .1);
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
            ylabel('xcorr');
            xlabel('time from lick s');
            axis tight
            hold off;
            
            %% excorr over shuff excorr cdf distr
            % relative swr from last lick
            sf2 = subaxis(2,2,2);
            anEXCSh = cell2mat([a.excorrShuf]);
            anEXC = [a.excorr];
            grps = [zeros(numel(anEXCSh),1); ones(numel(anEXC),1)];
            violin({anEXCSh, anEXC},...
                'facecolor',[.5 .5 .5; .6 .6 1;],'edgecolor','none');
            hold on
            b = boxplot([anEXCSh anEXC]',grps, ...
                'PlotStyle', 'compact', 'Symbol', '.','Color', 'k');
            set(gca,'XTickLabel',{' '})
            xticks([1 2])
            xticklabels({'control', 'real'})
            set(gca, 'FontSize', 10)
            legend off
            hold on;
            [p,h,stats] = ranksum(anEXCSh, anEXC);
            title(sprintf('ranksum p%.05f', p),'fontname','arial','fontsize', 10);
            ylabel('excess corr')
            hold off;
            
            %% polar distr phase clustering, swrLickPhase, meanMRVmag
            sf3 = subaxis(2,2,3);
            swrLPh = cell2mat({a.swrLickPhase}');
            swrLPhShuf = cell2mat([a.swrLickPhaseShuf]');
            [Rp, z] = circ_rtest(swrLPh);
            %         [kuippval, k, K] = circ_kuipertest(swrLPh, swrLPhShuff, 1, 1);
            polarhistogram(swrLPhShuf, 32, 'Normalization', 'pdf', 'edgealpha', .1,...
                'facecolor', [.5 .5 .5]);
            hold on
            polarhistogram(swrLPh, 32, 'Normalization', 'pdf', 'edgealpha', .1,...
                'facealpha', .8, 'facecolor', [.6 .6 1]);
            title(sprintf('rtest p%.04f', Rp))
            hold off
            axis tight
            
            %% phase mod
            sf4 = subaxis(2,2,4);
            anPHMODSh = cell2mat([a.phasemodShuf]);
            anPHMOD = [a.phasemod];
            grps = [zeros(numel(anPHMODSh),1); ones(numel(anPHMOD),1)];
            violin({anPHMODSh, anPHMOD},...
                'facecolor',[.5 .5 .5; .6 .6 1;],'edgecolor','none');
            hold on
            b = boxplot([anPHMODSh anPHMOD]',grps, ...
                'PlotStyle', 'compact', 'Symbol', '.','Color', 'k');
            set(gca,'XTickLabel',{' '})
            xticks([1 2])
            xticklabels({'control', 'real'})
            set(gca, 'FontSize', 10)
            legend off
            hold on;    
            [p,h,stats] = ranksum(anPHMODSh, anPHMOD);
            title(sprintf('ranksum p%.05f', p),'fontname','arial','fontsize', 10);
            ylabel('phasemod log(Z)')
            yl = ylim;
            ylim([yl(1) yl(2)+2])
            hold off;
            
            %% ----
            sprtit = sprintf('%s %s', animal, Fp.filtfunction(5:end));
            setSuperAxTitle(sprtit)
            if pausefigs; pause; end
            outdir = sprintf('%s/%s/', pconf.andef{4},Fp.filtfunction,animal);
            if savefigs; save_figure(outdir, sprtit, 'savefigas', savefigas); end
        end
    end
end