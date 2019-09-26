pconf = paramconfig;
create_filter = 0;
run_ff = 0;
save_ffdata = 0;
load_ffdata = 0;
combineEpochs = 0;
plotfigs = 1;
pausefigs = 0;
savefigs = 1;
savefigas = {'png'};
%% data filter params
Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_lickswrcorr';
Fp.params = {'savefigs', 'wtrackdays', Fp.filtfunction};
%% FF
Fp = load_filter_params(Fp);
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter,  ...
        'excludetime', Fp.timefilter,'iterator',Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff; F = runfilter(F);
    for d = 1:length(F); F(d).datafilter_params = Fp; end
end
if save_ffdata
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', sprintf('_%s', Fp.epochEnvironment));
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', sprintf('_%s', Fp.epochEnvironment));
end
% now that i have a few epoch done for D10, i was to bring them together..
% since i will now have several epochs per animal, and a measure for each epoch
% -- xcorr, excorr, mrvmag, phasemod
% i need to either show that all the animals are independently significant,
% or i need to combine them all, maybe individual epochs, into one pile and
% test that distribution vs the combination of all the shuffle
% distributions.. consult loren on this?
% really the best view would be to have
% - datachunks to check all the lick bouts and swr-associated lick bouts
% - examples from each animal (next)
% - each epoch  ( i have)
% - each animal (next)
% - all animals (next)
% - all animals different conditions?
%%
if combineEpochs
    smxc = {}; excorr = {}; mrvmag = {}; phmod = {}; smxcShuff = {}; excorrShuff = {};
    mrvmagShuff = {}; phmodShuff = {}; vecang = {};
    for ani = 1:length(F)
        % gather the real swr-lick measures across epochs
        smxc{ani} = cell2mat(cellfun(@(x) x.smthxc, F(ani).output, 'un', 0)');
        excorr{ani} = cell2mat(cellfun(@(x) x.excorr, F(ani).output, 'un', 0)');
        mrvmag{ani} = cell2mat(cellfun(@(x) x.meanMRVmag, F(ani).output, 'un', 0)');
        vecang{ani} = cell2mat(cellfun(@(x) x.vecang, F(ani).output, 'un', 0)');
        phmod{ani} = cell2mat(cellfun(@(x) x.phasemod, F(ani).output, 'un', 0)');
        swrLickPhase{ani} = cell2mat(cellfun(@(x) x.swrLickPhase, F(ani).output, 'un', 0)');
        % gather the swr-lick swr-lick shuffle distributions across epochs (1k per)
        smxcShuff{ani} = cell2mat(cellfun(@(x) x.smthxcShuff, F(ani).output, 'un', 0)');
        excorrShuff{ani} = cell2mat(cellfun(@(x) x.excorrShuff, F(ani).output, 'un', 0)');
        mrvmagShuff{ani} = cell2mat(cellfun(@(x) x.meanMRVmagShuff, F(ani).output, 'un', 0)');
%         vecangShuff{ani} = cellfun(@(x) x.vecangShuff, F(ani).output, 'un', 0)';
        phmodShuff{ani} = cell2mat(cellfun(@(x) x.phasemodShuff, F(ani).output, 'un', 0)');
        swrLickPhaseShuff{ani} = cell2mat(cellfun(@(x) x.swrLickPhaseShuff, F(ani).output, 'un', 0)');
    end
end
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
        time = cell2mat(cellfun(@(x) x.time, F(ani).output, 'un', 0)');
        try
        time = time(1,:);
        catch
            continue
        end
%         smxcshmean = nanmean(cell2mat(smxcShuff{ani}));
%         smxcshstd = std(cell2mat(smxcShuff{ani}));
%         plot(time, smxcshmean, 'color', [0 0 1 .2], 'linewidth', 1);
%         hold on;
%         fill([time'; flipud(time')],[smxcshmean'-smxcshstd';flipud(smxcshmean'+smxcshstd')],'b',...
%             'linestyle','none', 'facealpha', .1);
        % xcorr norm smooth
        smxcmean = nanmean(smxc{ani},1);
        smxcstd = nanstd(smxc{ani},[],1);
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
        anEXCSh = excorrShuff{ani};
        anEXCSh = anEXCSh(:);
        anEXC = excorr{ani};
        grps = [zeros(numel(anEXCSh),1); ones(numel(anEXC),1)];
        violin({anEXCSh, anEXC},...
            'facecolor',[.5 .5 .5; .6 .6 1;],'edgecolor','none'); 
        hold on
        b = boxplot([anEXCSh; anEXC],grps, ...
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
        
        swrLPh = swrLickPhase{ani}';
        swrLPhShuff = swrLickPhaseShuff{ani}';
%         anMRVmag = cell2mat(mrvmag{ani});
%         anvecang = cell2mat(vecang{ani});
%         anMRVmagShuff = cell2mat(mrvmagShuff{ani});
%         anvecangShuff = cell2mat(vecangShuff{ani});
        
        [Rp, z] = circ_rtest(swrLPh);
%         [kuippval, k, K] = circ_kuipertest(swrLPh, swrLPhShuff, 1, 1);
        polarhistogram(swrLPhShuff, 32, 'Normalization', 'pdf', 'edgealpha', .1,...
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
        anPHMODSh = phmodShuff{ani};
        anPHMODSh = anPHMODSh(:);
        anPHMOD = phmod{ani};
        grps = [zeros(numel(anPHMODSh),1); ones(numel(anPHMOD),1)];
        violin({anPHMODSh, anPHMOD},...
            'facecolor',[.5 .5 .5; .6 .6 1;],'edgecolor','none'); 
        hold on
        b = boxplot([anPHMODSh; anPHMOD],grps, ...
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
%         
%         histogram(phasemodShuff, 100, 'Normalization', 'probability', 'edgealpha', 0, 'facecolor', 'k');
%         hold on;
%         mrvsort = sort(phasemodShuff);
%         idxsig = round(sigperc*length(phasemodShuff));
%         line([mrvsort(idxsig) mrvsort(idxsig)], ylim, 'color', [0 0 0 .8], 'linestyle', '--');
%         hold on
%         line([phasemod phasemod], ylim, 'color', 'r');
%         modp = 1-sum(phasemod>phasemodShuff)/length(phasemodShuff);
%         title(sprintf('phasemod %.03f p%.03f Rpval%.03f', phasemod, modp, pval));
%         ylabel('probability')
%         xlabel('log(Rayleigh Z)')
%         axis tight
%         hold off
        

        %% ---- super axis -----
        sprtit = sprintf('%s lickswrcorr', animal);
        sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', gcf);
        iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
        set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment','center', ...
            'FontSize',12);
        h = get(gcf,'Children');
        set(gcf,'Children',flip(h)); % put super axis at bottom of axis stack. allows for zoom
        % ---- pause, save figs ----
        if pausefigs; pause; end
        outdir = sprintf('%s/lickswrcorr/', pconf.andef{4});
        if savefigs; save_figure(outdir, sprtit, 'savefigas', savefigas); end
    end
end
% animal
% plot a version of the perepoch fig, but for all the animals combined
% now that i have a distribution for the real, instead of just 1 number, i
% can do a boxwhisker between the real and the combined shuffled for each
% measure
