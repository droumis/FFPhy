
%{
remake final swr-lick xcorr figure for all animals.

remake reactivtion figure for all animals


%}

pconf = paramconfig;
create_filter = 1;
run_ff = 1;
load_ffdata = 0;
% combineEpochs = 0;

plotfigs = 1;
pausefigs = 0;
savefigs = 1;
savefigas = {'png', 'eps'};

%% data filter params
Fp.animals = {'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_lickswrcorr';
Fp.params = {'wtrackdays', Fp.filtfunction};

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
        'filetail', sprintf('_%s', Fp.epochEnvironment));
end

% %%
% if combineEpochs
%     smxc{ani} = 
%     
%     smxc = {}; excorr = {}; mrvmag = {}; phmod = {}; smxcShuff = {}; excorrShuff = {};
%     mrvmagShuff = {}; phmodShuff = {}; vecang = {};
%     for ani = 1:length(F)
%         gather the real swr-lick measures across epochs
%         smxc{ani} = cell2mat(cellfun(@(x) x.smthxc, F(ani).output, 'un', 0)');
%         excorr{ani} = cell2mat(cellfun(@(x) x.excorr, F(ani).output, 'un', 0)');
%         mrvmag{ani} = cell2mat(cellfun(@(x) x.meanMRVmag, F(ani).output, 'un', 0)');
%         vecang{ani} = cell2mat(cellfun(@(x) x.vecang, F(ani).output, 'un', 0)');
%         phmod{ani} = cell2mat(cellfun(@(x) x.phasemod, F(ani).output, 'un', 0)');
%         swrLickPhase{ani} = cell2mat(cellfun(@(x) x.swrLickPhase, F(ani).output, 'un', 0)');
%         gather the swr-lick swr-lick shuffle distributions across epochs (1k per)
%         smxcShuff{ani} = cell2mat(cellfun(@(x) x.smthxcShuff, F(ani).output, 'un', 0)');
%         excorrShuff{ani} = cell2mat(cellfun(@(x) x.excorrShuff, F(ani).output, 'un', 0)');
%         mrvmagShuff{ani} = cell2mat(cellfun(@(x) x.meanMRVmagShuff, F(ani).output, 'un', 0)');
%         vecangShuff{ani} = cellfun(@(x) x.vecangShuff, F(ani).output, 'un', 0)';
%         phmodShuff{ani} = cell2mat(cellfun(@(x) x.phasemodShuff, F(ani).output, 'un', 0)');
%         swrLickPhaseShuff{ani} = cell2mat(cellfun(@(x) x.swrLickPhaseShuff, F(ani).output, 'un', 0)');
%     end
% end

% plot a version of the perepoch fig, but for all epochs combined for each
%% plot

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
        time = cell2mat({a(1).time}');
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
            %% mag distr
        %     sf3 = subaxis(3,2,3);
        %     histogram(meanMRVmagShuff, 60, 'Normalization', 'probability', 'edgealpha', 0, 'facecolor', 'k');
        %     hold on;
        %     mrvsort = sort(meanMRVmagShuff);
        %     idxsig = round(sigperc*length(meanMRVmagShuff));
        %     line([mrvsort(idxsig) mrvsort(idxsig)], ylim, 'color', [0 0 0 .8], 'linestyle', '--');
        %     hold on
        %     line([meanMRVmag meanMRVmag], ylim, 'color', 'r');
        %     mrvp = 1-sum(meanMRVmag>meanMRVmagShuff)/length(meanMRVmagShuff);
        %     title(sprintf('MRVmag %.03f p%.03f', meanMRVmag, mrvp));
        %     ylabel('probability')
        %     axis tight
        %     hold off
%% polar distr phase clustering, swrLickPhase, meanMRVmag
        sf3 = subaxis(2,2,3);
        
        swrLPh = cell2mat({a.swrLickPhase}');
        swrLPhShuf = cell2mat([a.swrLickPhaseShuf]');
%         anMRVmag = cell2mat(mrvmag{ani});
%         anvecang = cell2mat(vecang{ani});
%         anMRVmagShuff = cell2mat(mrvmagShuff{ani});
%         anvecangShuff = cell2mat(vecangShuff{ani});
        
        [Rp, z] = circ_rtest(swrLPh);
%         [kuippval, k, K] = circ_kuipertest(swrLPh, swrLPhShuff, 1, 1);
        polarhistogram(swrLPhShuf, 32, 'Normalization', 'pdf', 'edgealpha', .1,...
            'facecolor', [.5 .5 .5]);
        hold on
        polarhistogram(swrLPh, 32, 'Normalization', 'pdf', 'edgealpha', .1,...
            'facealpha', .8, 'facecolor', [.6 .6 1]);
    	
%         polarplot([zeros(numel(swrLPh),1) swrLPh(:)]', ...
%             repmat([0 1],numel(swrLPh),1)', 'color', [0 0 0 .4], 'linewidth', 4);
%         a.ThetaAxisUnits = 'radians';
%         polarplot(anvecang, anMRVmag, '.', 'color', [0 0 0])
%         grid on
%         rticks([])
%         thetaticks([])
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
% animal
% plot a version of the perepoch fig, but for all the animals combined
% now that i have a distribution for the real, instead of just 1 number, i
% can do a boxwhisker between the real and the combined shuffled for each
% measure