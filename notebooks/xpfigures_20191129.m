
%{
continuation of eventtrig_20191117.. that was more about swr and lick time
based mod... this is focused on phase based mod. (but they need to be
combined?)

- beer is dmat maker compatible with output of forest.bear.cactus as well as barn.rat
- event trig LFP spect power: forest.bear.cactus.mushroom.leaf

- barn:rat is new, general Event-SU, with dmat:beer
    - barn:pig is old XP-SU time and phase mod
    - dfa_riptrigspiking is old SWR-SU

- barn:rat:beer:wheelbarrow is Event-SU timeMod (see eventtrig_20191117)
- barn:rat:beer:saw         is Event-SU phaseMod

- get spike containing ILI number within burst.. 
    - i had been within pig.. rat+saw took the place of pig.
    - also check that rat+saw gives the same result as pig

- space.alien is XP-SWR Time and Phase mod
- I need to run through space.alien and check its validity
    - i don't think i should squish XP-SWR into the barn.rat.saw/wheelbarrow pipeline..
        - because swr are pretty different than spikes
        - they have different iterators (space vs barn)
    - so for now I think keep space.alien, but maybe update it to make sure it's the 
        - phasemod shuffle phase/pct like in pig/saw (shuffle 0:100%)
        - e/xcorr shuffle time like in pig/wheelbarrow (shuffle -.5:.5s)

if im going to have conditions on swr-xp mod via, df alien.. i think i'll
need to 
either pass dmat 
%}

pconf = paramconfig;
eventTrigSWR = 1; % NOW.. PIPE:space.alien == per condition SWR mod
% eventTrigSpiking = 0; % PIPE:barn.rat.beer.saw/wheelbarrow == per condition/area Spike mod
% eventTrigLFP = 0; % PIPE:forest.bear.cactus.mushroom.beer.leaf == per condition/area LFP mod

eventType = 'lick'; %lick swr
% run FF
create_filter = 0;
run_ff = 0;
load_ffdata = 1;

%% DESIGN MAT MAKER.. WORKS WITH SPIKE AND LFP
% make_expvarCat = 0; % beer
% load_expvarCat = 0;

% %% Spike Event PhaseModulation
% calcSUphasemod = 0;
% loadSUPhaseMod = 0;
% % Gather across animals, per area, condition
% gatherPhaseModResults = calcSUphasemod;
gather_XP_mod_SWR = 1;
%% plot
plotfigs = 1;
showfigs = 1;
pausefigs = 0;
savefigs = 1;
savefigas = {'pdf'};

% XP phmod SWR
plot_XP_mod_SWR_perAn = 1;
plot_XP_mod_SWR_AllAn = 1;

% XP phmod SU
% plotSUphasemod = 0; % TODO?
% plotSpikePhaseModHeatRast = 0;
% plotCdfPolar = 0;

%% FF Data
Fp = [];
Fp.animals = {'JZ4'}; %'D10', 'D12', 'D13', 'JZ1', 'JZ3', 'JZ4'};
Fp.areas = {{'ca1', 'd'}, {'mec', 'deep'}, {'mec', 'supf'}};

% if eventTrigLFP
%     Fp.filtfunction = 'dfa_eventTrigLFP'; % Bear
%     if strcmp(eventType, 'lick')
%         expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
%         Fp.Label = 'wtrackLickTrigLFP';
%         Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
%             'excludeAfterLastWell', 'referenced', '4-350Hz',  'lickboutlicks', Fp.Label, Fp.filtfunction};
%     elseif strcmp(eventType, 'swr')
%         expvars = {'all', 'lickbouts', 'nolickbouts'};
%         Fp.Label = 'wtrackSWRTrigLFP';
%         Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
%             'excludeAfterLastWell', 'referenced', '4-350Hz',  'ripples', ...
%             Fp.Label, Fp.filtfunction};
%     end
%     
% elseif eventTrigSpiking
%     Fp.filtfunction = 'dfa_eventTrigSpiking'; % Rat
%     if strcmp(eventType, 'lick')
%         expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
%         Fp.Label = 'wtrackLickTrigSpiking';
%         Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
%             'excludeAfterLastWell', 'nonMU_cells', 'lickboutlicks', Fp.Label, Fp.filtfunction};
%     elseif strcmp(eventType, 'swr') %'excludeNoise',
%         expvars = {'all', 'lickbouts', 'nolickbouts'};
%         Fp.Label = 'wtrackSWRTrigSpiking';
%         Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
%             'excludeAfterLastWell', 'nonMU_cells', 'ripples', ...
%             Fp.Label, Fp.filtfunction}; % 'excludeNoise',
%     end
if eventTrigSWR
    Fp.filtfunction = 'dfa_lickswrcorr'; % city.alien
    if strcmp(eventType, 'lick')
%         expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
        Fp.Label = 'wXPTrigSWR';
        Fp.params = {'wtrackdays', 'excludePriorFirstWell', 'excludeAfterLastWell', ...
            'ripples>2', Fp.Label, Fp.filtfunction};
    else
        error('not yet implemented. eventTrigSWR is currently specific to XPtrigSWR')
    end
end
Fp = load_filter_params(Fp);
% if eventTrigLFP
%     wp = getWaveParams(Fp.waveSet);
% end
%%
if create_filter
%     if eventTrigLFP
%         F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', ...
%             Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
%     elseif eventTrigSpiking
%         F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes',...
%             Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator, 'cells',...
%             Fp.cellfilter);
    if eventTrigSWR
        F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'excludetime', ...
            Fp.timefilter, 'iterator', Fp.iterator);
    end
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F), 'un', 1);
    F = runfilter(F);
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.Label]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', ['_' Fp.Label]);
end
% %% FOR SPIKES make design mat to slice trials
% if make_expvarCat
%     if eventTrigSpiking
%         data = [];
%         for a = 1:length(F)
%             data(a).animal = F(a).animal;
%             idx = cell2mat({F(a).output{1}.index}');
%             [~, dayUnqIdx] = unique(idx(:,1));
%             days = idx(dayUnqIdx,[1]);
%             eps = idx(dayUnqIdx,[4:5]);
%             numEvPerDay = cell2mat({F(a).output{1}(dayUnqIdx).numEventsPerEp})';
%             data(a).day = [];
%             data(a).epoch = [];
%             for d = 1:size(eps,1)
%                 for e = 1:size(eps,2)
%                     data(a).day = [data(a).day; repmat(days(d), numEvPerDay(d,e), 1)];
%                     data(a).epoch = [data(a).epoch; repmat(eps(d,e), numEvPerDay(d,e),1)];
%                 end
%             end
%             data(a).evStart = cell2mat({F(a).output{1}(dayUnqIdx).eventTimes}');
%         end
%     else
%         data = lfpstack;
%     end
%     dmat = makeExpvarCatDesignMat(data, expvars, 'eventType', Fp.eventType);
% end
% 
% if load_expvarCat
%     outdir = 'expvarCat';
%     outpath = [pconf.andef{2},outdir,'/'];
%     dmat = load_data(outpath, [outdir,'_',Fp.env,'_',Fp.eventType], Fp.animals);
% end
% 
% %% Spike Phase mod
% if calcSUphasemod
%     pmodF = calcPhaseMod(F, dmat); % Saw
%     save_data(pmodF, 'results', [Fp.Label '_phasemod']);
% end
% if loadSUPhaseMod
%     pmodF = load_data('results', [Fp.Label '_phasemod'], Fp.animals);
% end
% 
% %% SPIKE Gather all animals Phasemod su per area, eventSet
% if gatherPhaseModResults
%     areaPhasemod = {};
%     allAnPmodF = cell2mat(arrayfun(@(x) pmodF(x).output{1}', 1:length(pmodF), 'un', 0)');
%     numESet = length(pmodF(1).dmatIdx);
%     
%     for ar = 1:length(Fp.areas) % per area
%         % find cells in this area
%         areaIdx = strcmp({allAnPmodF.area}', Fp.areas{ar}{1});
%         subareaIdx = ~cellfun(@isempty, strfind({allAnPmodF.subarea}', Fp.areas{ar}{2}), 'un', 1);
%         iareaIdx = find(all([areaIdx subareaIdx],2));
%         iareaIdx = iareaIdx(arrayfun(@(x) ~isempty(allAnPmodF(x).phasemod), ...
%             iareaIdx,'un',1));
%         
%         for iv = 1:numESet % per eventSet
%             gud = [];
%             dumpy = [];
%             phasemod = [];
%             mPctChangeSh = [];
%             for i = 1:length(iareaIdx)
%                 try
%                     m = allAnPmodF(iareaIdx(i)).phasemod{iv};
%                     %                     mSh = allAnPmodF(iareaIdx(i)).mPctChangeSh{iv};
%                     if ~isempty(m)
%                         phasemod = [phasemod; m];
%                         mPctChangeSh = [mPctChangeSh; mSh];
%                         gud = [gud; i];
%                     end
%                 catch
%                     continue
%                 end
%             end
%             areaPhasemod{ar,iv} = phasemod;
%         end
%     end
% end
%% Gather the data from city.squirrel F needed for plotting
if gather_XP_mod_SWR
    NormxcMean = {};
    NormxcZmean = {};
    NormxcZsem = {};
    swrLPh = {};
    swrLPhShuf = {};
    anPHMODSh = {};
    anPHMOD = {};
    meanMRVmag = {};
    vecang = {};
    
    for a = 1:length(F)
        idata = [F(a).output{:}];
        nonempty = find(~cell2mat(cellfun(@isempty, {idata.time}, 'un', 0)));
        time = cell2mat({idata(nonempty(1)).time}');
        
        % xcorr norm
        nxc = cell2mat({idata.normxc}');
        NormxcMean{a,1} = nanmean(nxc, 1);
        %     NormxcSem{a,1} = sem(nxc,1);
        
        % z score smooth xc
        dNormxcZ = {};
        for d = 1:length(idata) % for each day
            dNormxc = idata(d).smthxc;
            dNormxcSh = idata(d).smthxcSh;
            dNormxcZ{d} = (dNormxc - nanmean(dNormxcSh,1)) ./ nanstd(dNormxcSh,[],1);
        end
        NormxcZ = cat(1,dNormxcZ{:});
        NormxcZsem{a,1} = sem(NormxcZ,1);
        NormxcZmean{a,1} = mean(NormxcZ,1);
        
        % Polar Clustering
        swrLPh{a,1} = cell2mat(cat(1,{idata.swrLickPhase}'));
        swrLPhShuf{a,1} = cell2mat(cat(1,{idata.swrLickphaseSh}'));
        [Rp, z] = circ_rtest(swrLPh{a});
%         cell2mat(cat(1,{idata.meanMRVmag}'));
        meanvec = mean(exp(1i*swrLPh{a})); % get mean resultant vector
        meanMRVmag{a,1} = abs(meanvec); % vector magnitude
        vecang{a,1} = angle(meanvec);
        
        % phasemod stingray
        for i = 1:size(swrLPhShuf{a},2)
            [~, z] = circ_rtest(swrLPhShuf{a}(:,i));
            anPHMODSh{a,1}(i) = log(z);
        end
        anPHMOD{a,1} = [idata.phasemod];
    end
    
    NormxcMeanAll = nanmean(cat(1,NormxcMean{:}),1);
    NormxcSemAll = sem(cat(1,NormxcMean{:}),1);
    
    NormxcZmeanAll = nanmean(cell2mat(NormxcZmean),1);
    NormxcZsemAll = sem(cell2mat(NormxcZsem),1);
    
    swrLPhAll = cell2mat(swrLPh);
    swrLPhShufAll = cell2mat(swrLPhShuf);
    
    anPHMODAll = cell2mat(anPHMOD'); % cat all days from all ans
    anPHMODShAll = cell2mat(anPHMODSh');
end
%% PLOT=====================================================================

if plotfigs
    % Plot Event1 mod Event2
    if plot_XP_mod_SWR_perAn
        figname = 'wXPmodSWR';
        Pp=load_plotting_params({'defaults','dfa_lickswrcorr', figname});
        for a = 1:length(F)
            animal = F(a).animal{3};
            ifig = init_plot(showfigs, Pp.position);
            %% Xcorr norm bar hist
            sf1 = subaxis(2,2,1,Pp.posparams{:});
            sf1.Tag = 'normxc';
            bar(time, NormxcMean{a},'FaceColor','k')
%             errorbar(time,NormxcMean,NormxcSem,'.','Color',[0.5 0.5 0.5])
            % Pp
            axis tight
            title(' MeanNorm perDay')
            xlabel('Time of SWRs relative to XPs')
            ylabel({'Norm xcorr'})
            xlim([-0.5 0.5])
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
            set(gca,'xtick',[-0.5:0.1:.5])
            hold off;
            
            %% smooth Xcorr z scored
            sf1 = subaxis(2,2,3,Pp.posparams{:});
            sf1.Tag = 'smthxcZ';
            plot(time, NormxcZmean{a}', 'k')
            hold on
            fill([time'; flipud(time')],[NormxcZmean{a}'-NormxcZsem{a}'; ...
                flipud(NormxcZmean{a}'+NormxcZsem{a}')],'k',...
                'linestyle','none', 'facealpha', .2);

            axis tight
            title(' MeanZscoredSmth perDay relShuf')
            xlabel('Time of SWRs relative to XPs')
            ylabel({'Z-scored xcorr'})
            xlim([-0.5 0.5])
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
            set(gca,'xtick',[-0.5:0.1:.5])
            hold off;
            %%
%             %% excorr over shuff excorr cdf distr
%             % relative swr from last lick
%             sf2 = subaxis(2,2,3);
%             anEXCSh = [idata.excorrSh];
%             anEXC = [idata.excorr];
%             grps = [zeros(numel(anEXCSh),1); ones(numel(anEXC),1)];
%             violin({anEXCSh, anEXC},...
%                 'facecolor',[.5 .5 .5; .6 .6 1;],'edgecolor','none');
%             hold on
%             b = boxplot([anEXCSh anEXC]',grps, ...
%                 'PlotStyle', 'compact', 'Symbol', '.','Color', 'k');
%             set(gca,'XTickLabel',{' '})
%             xticks([1 2])
%             xticklabels({'control', 'real'})
%             set(gca, 'FontSize', 10)
%             legend off
%             hold on;
%             [p,h,stats] = ranksum(anEXCSh, anEXC);
%             title(sprintf('ranksum p%.05f', p),'fontname','arial','fontsize', 10);
%             ylabel('excess corr')
%             hold off;
%             
            %% polar distr phase clustering, swrLickPhase, meanMRVmag
            sf3 = subaxis(2,2,2);
%             swrLPh = cell2mat({idata.swrLickPhase}');
%             swrLPhShuf = cell2mat({idata.swrLickphaseSh}');
%             [Rp, z] = circ_rtest(swrLPh);
            %         [kuippval, k, K] = circ_kuipertest(swrLPh, swrLPhShuff, 1, 1);
            k = .7;
            polarhistogram(swrLPhShuf{a}(:), 32, 'Normalization', 'pdf', 'edgealpha', .1,...
                'facecolor', [k k k], 'DisplayName', 'shuffled');
            hold on
            polarhistogram(swrLPh{a}, 32, 'Normalization', 'pdf', 'edgealpha', .1,...
                'facealpha', .9, 'facecolor', [.5 .5 .9], 'DisplayName', 'phasemod /Day');
            legend show
            
            thetaticks(0:45:315)
            pax = gca;
            pax.ThetaAxisUnits = 'radians';
            pax.ThetaZeroLocation = 'left';
            b = .3;
            pax.RColor = [b b b];
            pax.ThetaColor = [b b b];
            title(sprintf('rtest p%.04f', Rp))
            hold off
            axis tight
            
            %% phase mod
            sf4 = subaxis(2,2,4);
            
%             for i = 1:size(swrLPhShuf,2)
%               [~, z] = circ_rtest(swrLPhShuf(:,i));
%               anPHMODSh(i) = log(z);
%             end
%             anPHMODSh = cell2mat([idata.swrLickphaseSh]);
%             anPHMOD = [idata.phasemod];
            
            grps = [zeros(numel(anPHMODSh{a}),1); ones(numel(anPHMOD{a}),1)];
            violin({anPHMODSh{a}, anPHMOD{a}},...
                'facecolor',[k k k; .5 .5 .9;],'edgecolor','none');
            hold on
            b = boxplot([anPHMODSh{a} anPHMOD{a}]',grps, ...
                'PlotStyle', 'compact', 'Symbol', '.','Color', [.2 .2 .2]);
            set(gca,'XTickLabel',{' '})
            xticks([1 2])
            xticklabels({'shuffled', 'phasemod /Day'})
            set(gca, 'FontSize', 10)
            legend off
            hold on;
            [p,h,stats] = ranksum(anPHMODSh{a}, anPHMOD{a});
            title(sprintf('ranksum p%.05f', p),'fontname','arial','fontsize', 10);
            ylabel('phasemod log(Z)')
            yl = ylim;
            ylim([yl(1) yl(2)+2])
            hold off;
            
            %% ----
            stit = sprintf('%s %s', figname, animal);
            setSuperAxTitle(stit);
            if pausefigs
                pause;
            end
            if savefigs
                save_figure([pconf.andef{4} '/' figname], stit, 'savefigas', ...
                    savefigas);
            end
        end
    end
    
        %% Plot All An. Event1 mod Event2
    if plot_XP_mod_SWR_AllAn
        figname = 'wXPmodSWRall';
        Pp=load_plotting_params({'defaults', figname});
        ifig = init_plot(showfigs, Pp.position);
        
        %% Xcorr norm
        sf1 = subaxis(2,2,1,Pp.posparams{:});
        bar(time, NormxcMeanAll,'FaceColor','k')
        % Pp
        axis tight
        title(' MeanNorm perDay')
        xlabel('Time of SWRs relative to XPs')
        ylabel({'Norm xcorr'})
        xlim([-0.5 0.5])
        line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
        set(gca,'xtick',[-0.5:0.1:.5])
        hold off;
        
        %% smooth Xcorr z scored
        sf1 = subaxis(2,2,3,Pp.posparams{:});
        plot(time, NormxcZmeanAll', 'k')
        hold on
        fill([time'; flipud(time')],[NormxcZmeanAll'-NormxcZsemAll'; ...
            flipud(NormxcZmeanAll'+NormxcZsemAll')],'k',...
            'linestyle','none', 'facealpha', .2);
        
        axis tight
        title(' MeanZscoredSmth perDay relShuf')
        xlabel('Time of SWRs relative to XPs')
        ylabel({'Z-scored xcorr'})
        xlim([-0.5 0.5])
        line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
        set(gca,'xtick',[-0.5:0.1:.5])
        hold off;
        
        %% polar distr phase clustering, swrLickPhase, meanMRVmag
        sf3 = subaxis(2,2,2);
        k = .7;
        polarhistogram(swrLPhShufAll(:), 32, 'Normalization', 'pdf', 'edgealpha', .1,...
            'facecolor', [k k k], 'DisplayName', 'shuffled');
        hold on
        polarhistogram(swrLPhAll, 32, 'Normalization', 'pdf', 'edgealpha', .1,...
            'facealpha', .9, 'facecolor', [.5 .5 .9], 'DisplayName', 'phasemod /Day');
        legend show
        polarscatter(cell2mat(vecang), cell2mat(meanMRVmag), 'k.');
        
        thetaticks(0:45:315)
        pax = gca;
        pax.ThetaAxisUnits = 'radians';
        pax.ThetaZeroLocation = 'left';
        b = .3;
        pax.RColor = [b b b];
        pax.ThetaColor = [b b b];
        title(sprintf('rtest p%.04f', Rp))
        hold off
        axis tight
        
        %% phase mod
        sf4 = subaxis(2,2,4);
        
        grps = [zeros(numel(anPHMODShAll),1); ones(numel(anPHMODAll),1)];
        violin({anPHMODShAll, anPHMODAll},...
            'facecolor',[k k k; .5 .5 .9;],'edgecolor','none');
        hold on
        b = boxplot([anPHMODShAll anPHMODAll]',grps, ...
            'PlotStyle', 'compact', 'Symbol', '.','Color', [.2 .2 .2]);
        set(gca,'XTickLabel',{' '})
        xticks([1 2])
        xticklabels({'shuffled', 'phasemod /Day'})
        set(gca, 'FontSize', 10)
        legend off
        hold on;
        [p,h,stats] = ranksum(anPHMODShAll, anPHMODAll);
        title(sprintf('ranksum p%.05f', p),'fontname','arial','fontsize', 10);
        ylabel('phasemod log(Z)')
        yl = ylim;
        ylim([yl(1) yl(2)+2])
        hold off;
        
        %% ----
        stit = sprintf('%s AllAn', figname);
        setSuperAxTitle(stit);
        if pausefigs
            pause;
        end
        if savefigs
            save_figure([pconf.andef{4} '/' figname], stit, 'savefigas', ...
                savefigas);
        end
    end
end
    %% plot examples
    % for every lick bout, plot SU, group by area
    
    %% plot phasemod per unit
    % time axis and polar axis?
    % polar axis
    
    % 1 spike per line.. this will get crowded.. 
    % maybe use a histogram instead
    
%     %% plot XP-phase of heatrast per animal, then all animal
%     if plotSpikePhaseModHeatRast
%         figname = 'suLickPhaseModHeatRastWtrack';
%         Pp=load_plotting_params({'defaults', figname});
%         hbins = round(Pp.numBins/2);
%         binphase = linspace(0, 2*pi, Pp.numBins);
%         binphC = binphase(1:end-1) + diff(1:2)/2; % bin centers
%         allAnSpikeXPhHist = {};
%         for a = 1:length(pmodF)
%             animal = pmodF(a).animal{3};
%             fprintf('%s\n', animal);
%             for ar = 1:length(Fp.areas)
%                 ifig = init_plot(showfigs, Pp.position);
%                 % find cells in this area
%                 areaIdx = strcmp({pmodF(a).output{1}.area}', Fp.areas{ar}{1});
%                 subareaIdx = ~cellfun(@isempty, strfind({pmodF(a).output{1}.subarea}', Fp.areas{ar}{2}));
%                 iareaIdx = find(all([areaIdx subareaIdx],2));
%                 for iv = 1:length(pmodF(a).dmatIdx)
%                     %% all su in area heatraster  
%                     sf = subaxis(1,length(pmodF(a).dmatIdx),iv, Pp.posparams{:});
%                     sf.Tag = 'heatrast';
%                     % make heatraster from all the clusters in this area, condition
%                     spikeXPhHist = {};
%                     for f = 1:length(iareaIdx)
%                         try
%                         spikeXPhHist{f,1} = histcounts(pmodF(a).output{1}(iareaIdx(f)...
%                             ).spikeLickPhase{iv}, binphase); %,'Normalization', 'probability');
%                         catch
%                             continue
%                         end
%                     end
%                     spikeXPhHist = cell2mat(spikeXPhHist);
%                     % pad edges
%                     padS = spikeXPhHist(:,hbins+1:end);
%                     padE = spikeXPhHist(:,1:hbins);
%                     spkXPHistPad = [padS spikeXPhHist padE];
%                     % sort based on sm max idx
%                     spkXPHistPadSm = smoothdata(spkXPHistPad, 2, 'loess', round(size(spkXPHistPad,2)/2));
%                     [~, imx] = max(spkXPHistPadSm(:,hbins:end-hbins),[],2);
%                     [~, modsortIdx] = sort(imx, 1, 'descend');
%                     spkXPHistPadSmSort = spkXPHistPadSm(modsortIdx,:);
%                     spkXPHistPadSmSortZ = zscore(spkXPHistPadSmSort, [],2);
%                     % plot
%                     imagesc(padBinPhC, [1:size(spkXPHistPadSmSortZ,1)], spkXPHistPadSmSortZ);
%                     xlim([-1*pi 3*pi])
%                     xticks([-1*pi 0 pi 2*pi 3*pi])
%                     xticklabels({'-1\pi', '0','\pi','2\pi','3\pi'})
%                     caxis(sf,[-1.25 2.25])
%                     colormap(magma);
%                     hcb = colorbar;
%                     hax = gca;
%                     hax.YTickLabel = flipud(hax.YTickLabel);
%                     set(get(hcb,'ylabel'),'String','z-score firing rate', 'rotation', 270, 'fontsize', 12)
%                     ylabel(sprintf('%s SU num', strjoin(Fp.areas{:,ar})))
%                     title(sprintf('%s', pmodF(a).dmatIdx{iv}))
%                     xlabel('lick phase')
%                     allAnSpikeXPhHist{ar, iv, a} = spikeXPhHist;
%                 end
%                 stit = sprintf('%s %s %s', figname, animal, strjoin(Fp.areas{:,ar}));
%                 setSuperAxTitle(stit);
%                 if pausefigs
%                     pause;
%                 end
%                 if savefigs
%                     save_figure([pconf.andef{4} '/' figname '/' animal], stit, 'savefigas', ...
%                         savefigas);
%                 end
%             end
%         end
%         %% all an, per area, condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         for ar = 1:size(allAnSpikeXPhHist,1)
%             ifig = init_plot(showfigs, Pp.position);
%             for iv = 1:size(allAnSpikeXPhHist,2)
%                 %%
%                 sf = subaxis(1,size(allAnSpikeXPhHist,2),iv, Pp.posparams{:});
%                 sf.Tag = 'heatrast';
%                 spikeXPhHist = cell2mat(squeeze(allAnSpikeXPhHist(ar,iv,:)));
%                 
%                 % pad edges
%                 padS = spikeXPhHist(:,hbins+1:end);
%                 padE = spikeXPhHist(:,1:hbins);
%                 spkXPHistPad = [padS spikeXPhHist padE];
%                 % sort based on sm max idx
%                 spkXPHistPadSm = smoothdata(spkXPHistPad, 2, 'loess', round(size(spkXPHistPad,2)/2.7));
%                 [~, imx] = max(spkXPHistPadSm(:,hbins:end-hbins),[],2);
%                 [~, modsortIdx] = sort(imx, 1, 'descend');
%                 spkXPHistPadSmSort = spkXPHistPadSm(modsortIdx,:);
%                 spkXPHistPadSmSortZ = zscore(spkXPHistPadSmSort, [],2);
%                 % plot
%                 imagesc(padBinPhC, [1:size(spikeXPhHist,1)], spkXPHistPadSmSortZ);
%                 xlim([-1*pi 3*pi])
%                 xticks([-1*pi 0 pi 2*pi 3*pi])
%                 xticklabels({'-1\pi', '0','\pi','2\pi','3\pi'})
%                 caxis(sf,[-1.1 2.1])
%                 colormap(magma);
%                 hcb = colorbar;
%                 hax = gca;
%                 hax.YTickLabel = flipud(hax.YTickLabel);
%                 set(get(hcb,'ylabel'),'String','z-score firing rate', 'rotation', 270, 'fontsize', 12)
%                 ylabel(sprintf('%s SU num', strjoin(Fp.areas{:,ar})))
%                 title(sprintf('%s', pmodF(a).dmatIdx{iv}))
%                 xlabel('lick phase')
%                 %%
%             end
%             stit = sprintf('%s AllAn %s', figname, strjoin(Fp.areas{:,ar}));
%             setSuperAxTitle(stit);
%             if pausefigs
%                 pause;
%             end
%             if savefigs
%                 save_figure([pconf.andef{4} '/' figname], stit, 'savefigas', savefigas);
%             end
%         end
%     end
%     
%     %% plot the mrv of the sig phasemod su %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if plotCdfPolar
%         figname = 'wXPSU-cdfPolar';
%         allAnPhasemod= {};
%         allAnVecang = {};
%         allAnPhasemodSh = {};
%         allAnVecangSh = {};
%         for a = 1:length(pmodF)
%             animal = pmodF(a).animal{3};
%             Pp = load_plotting_params({'defaults', figname});
%             for ar = 1:length(Fp.areas)
%                 ifig = init_plot(showfigs, Pp.position);
%                 % find SU from this Area
%                 areaIdx = ~cellfun(@isempty, strfind({pmodF(a).output{1}.area}', Fp.areas{ar}{1}));
%                 subareaIdx = ~cellfun(@isempty, strfind({pmodF(a).output{1}.subarea}', Fp.areas{ar}{2}));
%                 iareaIdx = find(all([areaIdx subareaIdx],2));
%                 for iv = 1:length(pmodF(a).dmatIdx)
%                     % collect data for this animal, area, condition
%                     phasemod = {};
%                     vecang = {};
% %                     mrvmag = {};
%                     phasemodSh = {};
%                     vecangSh = {};
% %                     mrvmagSh = {};                    
%                     for spikeXPhHist = 1:length(iareaIdx)
%                         try
%                             phasemod{spikeXPhHist} = pmodF(a).output{1}(iareaIdx(spikeXPhHist)).phasemod{iv};
%                             vecang{spikeXPhHist} = pmodF(a).output{1}(iareaIdx(spikeXPhHist)).vecang{iv};
% %                             mrvmag{spikeXPhHist} = pmodF(a).output{1}(iareaIdx(spikeXPhHist)).meanMRVmag{iv};
% %                             pval{spikeXPhHist} = pmodF(a).output{1}(iareaIdx(spikeXPhHist)).pval{iv};
%                             phasemodSh{spikeXPhHist} = pmodF(a).output{1}(iareaIdx(spikeXPhHist)).phasemodSh{iv};
%                             vecangSh{spikeXPhHist} = pmodF(a).output{1}(iareaIdx(spikeXPhHist)).vecangSh{iv};
% %                             mrvmagSh{spikeXPhHist} = pmodF(a).output{1}(iareaIdx(spikeXPhHist)).meanMRVmagSh{iv};
%                         catch
%                             continue
%                         end
%                     end
%                     phasemod = cell2mat(phasemod)';
%                     vecang = cell2mat(vecang)';
% %                     mrvmag = cell2mat(mrvmag)';
%                     phasemodSh = cell2mat(phasemodSh)';
%                     vecangSh = cell2mat(vecangSh)';
% %                     mrvmagSh = cell2mat(mrvmagSh)';
%                     allAnPhasemod{ar, iv, a} = phasemod;
%                     allAnVecang{ar, iv, a} = vecang;
%                     allAnPhasemodSh{ar, iv, a} = phasemodSh;
%                     allAnVecangSh{ar, iv, a} = vecangSh;
%                     %% cdf SU phasemod
%                     sf = subaxis(2,length(pmodF(a).dmatIdx),iv,Pp.posparams{:});
%                     sf.Tag = strjoin(Fp.areas{:,ar});
%                     [h, b] = histcounts(phasemodSh, Pp.logZBins, 'Normalization','cdf');
%                     logZticks = b(1:end-1)+diff(b)/2;
%                     plot(logZticks, h*100, 'k', 'DisplayName', 'shuffle')
%                     hold on
%                     [h, b] = histcounts(1-phasemod, Pp.logZBins, 'Normalization','cdf');
%                     logZticks = b(1:end-1)+diff(b)/2;
%                     plot(logZticks, h*100, 'DisplayName', 'CA1')
%                     hold on
%                     axis tight
%                     legend('Location','northwest','Orientation', 'vertical')
%                     
%                     phShSigLen = length(1-phasemodSh) * Pp.sigPct;
%                     phShufSort = sort(phasemodSh);
%                     sigThresh = phShufSort(round(phShSigLen));
%                     line([sigThresh sigThresh], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', 1, ...
%                         'HandleVisibility','off')
% %                                 ax = gca;
% %                                 ax.YDir = 'reverse'; % to match karalis
%                     xlabel('phase-modulation')
%                     ylabel('% units');
%                     grid on
%                     xl = xlim;
%                     xlim([-4.5 xl(2)])
%                     hold off
%                     title(sprintf('%s', pmodF(a).dmatIdx{iv}));
%                     %% sig phasemod mrv POLAR
%                     sf = subaxis(2,length(pmodF(a).dmatIdx),iv+length(pmodF(a).dmatIdx),Pp.posparams{:});
% %                     sf = subaxis(1,1,2,Pp.posparams{:});
%                     sf.Tag = strjoin(Fp.areas{:,ar});
% 
%                     idxSig = find(phasemod > sigThresh);
%                     idxNSig = find(phasemod < sigThresh);
% %                     [thSSh, ~] = sort(vecangSh);
% %                     [thN, ~] = sort(vecang(idxNSig));
% %                     [thS, thIdx] = sort(vecang(idxSig));
% %                     r = mrvmag(idxSig);
% %                     rS = r(thIdx);
% %                     rsAll = mrvmag(thIdxAll);
%                     polarhistogram(vecangSh, 24, 'Normalization', 'pdf','facecolor', [0 0 0])
%                     hold on
%                     polarhistogram(vecang(idxNSig), 24, 'Normalization', 'pdf', 'facecolor', [0 0 1], 'facealpha', .5); %, 'markerfacecolor', [0 0 0 .5])
%                     polarhistogram(vecang(idxSig), 24, 'Normalization', 'pdf', 'facecolor', [1 0 0], 'facealpha', .5); %, 'markerfacecolor', [0 0 0 .5])                    
% %                     hold on
% %                     polarscatter(thSAll, rsAll, 20, 'filled', 'markerfacealpha', .5, 'markerfacecolor', [0 0 0])
% %                     polarscatter(thS, rS, 20, 'filled', 'markerfacealpha', .5); %, 'markerfacecolor', [0 0 0 .5])
% %                     hold on
%                     hold off
%                     Ax = gca; % current axes
%                     %             Ax.ThetaGrid = 'off';
%                     Ax.RGrid = 'on';
%                     Ax.ThetaGrid = 'on';
%                     thetaticks([0 45 90 135 180 225 270 315])
%                     Ax.ThetaAxisUnits = 'radians';
%                     Ax.RTickLabel = [];
%                     Ax.ThetaAxis.Label.String = 'MRV sigSU';
%                     Ax.ThetaAxis.Label.Rotation = 270;
%                     axis tight
% %                     Ax.ThetaTickLabel = [];
%                 end
%                 %%
%                 stit = sprintf('%s %s %s %s', figname, animal, Fp.env, strjoin(Fp.areas{:,ar}));
%                 setSuperAxTitle(stit);
%                 if pausefigs
%                     pause
%                 end
%                 if savefigs
%                     strsave = save_figure([pconf.andef{4} '/' figname '/' animal], stit);
%                     strsave = save_figure([pconf.andef{4} '/' figname '/' animal], stit);
% %                     F(a).figs.spikePhaseCumPolar{ar} = strsave;
%                 end
%             end
%         end
%          %% all an, per area, condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         for ar = 1:size(allAnPhasemod,1)
%             ifig = init_plot(showfigs, Pp.position);
%             for iv = 1:size(allAnPhasemod,2)
%                 phasemod = cell2mat(squeeze(allAnPhasemod(ar,iv,:)));
%                 vecang = cell2mat(squeeze(allAnVecang(ar,iv,:)));
%                 phasemodSh = cell2mat(squeeze(allAnPhasemodSh(ar,iv,:)));
%                 vecangSh = cell2mat(squeeze(allAnVecangSh(ar,iv,:)));
%                 
%                 %% cdf SU phasemod
%                 sf = subaxis(2,length(pmodF(a).dmatIdx),iv,Pp.posparams{:});
%                 sf.Tag = strjoin(Fp.areas{:,ar});
%                 [h, b] = histcounts(phasemodSh, Pp.logZBins, 'Normalization','cdf');
%                 logZticks = b(1:end-1)+diff(b)/2;
%                 plot(logZticks, h*100, 'k', 'DisplayName', 'shuffle')
%                 hold on
%                 [h, b] = histcounts(1-phasemod, Pp.logZBins, 'Normalization','cdf');
%                 logZticks = b(1:end-1)+diff(b)/2;
%                 plot(logZticks, h*100, 'DisplayName', strjoin(Fp.areas{:,ar}))
%                 hold on
%                 axis tight
%                 legend('Location','northwest','Orientation', 'vertical')
%                 
%                 phShSigLen = length(1-phasemodSh) * Pp.sigPct;
%                 phShufSort = sort(phasemodSh);
%                 sigThresh = phShufSort(round(phShSigLen));
%                 line([sigThresh sigThresh], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', 1, ...
%                     'HandleVisibility','off')
%                 %                                 ax = gca;
%                 %                                 ax.YDir = 'reverse'; % to match karalis
%                 xlabel('phase-modulation')
%                 ylabel('% units');
%                 grid on
%                 xl = xlim;
%                 xlim([-4.5 xl(2)])
%                 hold off
%                 title(sprintf('%s', pmodF(a).dmatIdx{iv}));
%                 
%                 
%                 %% sig phasemod mrv POLAR
%                 sf = subaxis(2,length(pmodF(a).dmatIdx),iv+length(pmodF(a).dmatIdx),Pp.posparams{:});
%                 %                     sf = subaxis(1,1,2,Pp.posparams{:});
%                 sf.Tag = strjoin(Fp.areas{:,ar});
%                 idxSig = find(phasemod > sigThresh);
%                 idxNSig = find(phasemod < sigThresh);
%                 %                     [thSSh, ~] = sort(vecangSh);
%                 %                     [thN, ~] = sort(vecang(idxNSig));
%                 %                     [thS, thIdx] = sort(vecang(idxSig));
%                 %                     r = mrvmag(idxSig);
%                 %                     rS = r(thIdx);
%                 %                     rsAll = mrvmag(thIdxAll);
%                 polarhistogram(vecangSh, 24, 'Normalization', 'pdf','facecolor', [0 0 0])
%                 hold on
%                 polarhistogram(vecang(idxNSig), 24, 'Normalization', 'pdf', 'facecolor', [0 0 1], 'facealpha', .5); %, 'markerfacecolor', [0 0 0 .5])
%                 polarhistogram(vecang(idxSig), 24, 'Normalization', 'pdf', 'facecolor', [1 0 0], 'facealpha', .5); %, 'markerfacecolor', [0 0 0 .5])
%                 %                     hold on
%                 %                     polarscatter(thSAll, rsAll, 20, 'filled', 'markerfacealpha', .5, 'markerfacecolor', [0 0 0])
%                 %                     polarscatter(thS, rS, 20, 'filled', 'markerfacealpha', .5); %, 'markerfacecolor', [0 0 0 .5])
%                 %                     hold on
%                 hold off
%                 Ax = gca; % current axes
%                 %             Ax.ThetaGrid = 'off';
%                 Ax.RGrid = 'on';
%                 Ax.ThetaGrid = 'on';
%                 thetaticks([0 45 90 135 180 225 270 315])
%                 Ax.ThetaAxisUnits = 'radians';
%                 Ax.RTickLabel = [];
%                 Ax.ThetaAxis.Label.String = 'MRV sigSU';
%                 Ax.ThetaAxis.Label.Rotation = 270;
%                 axis tight
%                 %%
%             end
%             stit = sprintf('%s AllAn %s', figname, strjoin(Fp.areas{:,ar}));
%             setSuperAxTitle(stit);
%             if pausefigs
%                 pause;
%             end
%             if savefigs
%                 save_figure([pconf.andef{4} '/' figname], stit, 'savefigas', savefigas);
%             end
%         end
%     end
% end