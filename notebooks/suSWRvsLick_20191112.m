%{
- the central question
- the su lick-resp is more complicated and rhtyhmic

ca1. su with principal unit trough/peak (how did mari do this?)

run eventTrigSpiking using SWR as events.
run lickXCorrSpikes using lick as events and specifically looking at
interlick intervals.. so the phasic stuff.

A:: heatraster of all clust swr-trig su sorted by resp max

B:: heatraster of lick phase su, sorted by the sorting in A

linregress the swr-mod score vs phase-mod
1: none: swr-resp neurons in
2: pos: swr-resp neurons are also lick-resp
3:

extra: spatial firing rate map
%}

pconf = paramconfig;
runSWRTrigSpiking = 0;
runLickPhaseSpiking = 1;

create_filter = 0;
run_ff = 0;
load_ffdata = 0;

% via dfa_eventTrigSpiking, calcMod. per animal clusts mean psth, stat
calcMod = 0;
loadMod = 0;
plotHeatRaster = 0;
% plotStatMod = 0;

% flat phase
plotClustILPC = 0;
% via dfa_lickXCorrSpikes
plotILPTHeatRaster = 0;
plotCumPolar = 1;

showfigs = 0;
pausefigs = 0;
savefigs = 1;
animals = {'D10', 'D13', 'JZ1', 'JZ4'};

%% SU SWR Response
if runSWRTrigSpiking
    Fp = [];
    Fp.Label = 'swrTrigSUmod';
    Fp.animals = animals;
    Fp.filtfunction = 'dfa_eventTrigSpiking';
    Fp.params = {Fp.Label, 'ripples', 'wtrackdays', 'valid_ntrodes', Fp.filtfunction, 'nonMU_cells'};
    Fp = load_filter_params(Fp);
    Fp.areas = {{'ca1'}, {'mec'}};
    Fp.filetail = sprintf('_%s_%s', Fp.eventType, Fp.env);
    %% FF
    if create_filter
        F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes',...
            Fp.tetfilter, 'excludetime', Fp.timefilter,'iterator',Fp.iterator, 'cells',...
            Fp.cellfilter);
        F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    end
    if run_ff
        F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F), 'un', 1);
        F = runfilter(F);
        save_data(F, [pconf.andef{2} Fp.Label], Fp.filtfunction, ...
            'filetail', Fp.filetail);
    end
    if load_ffdata
        F = load_data([pconf.andef{2} Fp.Label], Fp.filtfunction, Fp.animals, ...
            'filetail', Fp.filetail);
    end
    if calcMod
        modF = calcSUmod(F);
        save_data(modF, 'results', Fp.Label, 'filetail', Fp.filetail);
    end
    if loadMod
        modF = load_data('results', Fp.Label, Fp.animals, 'filetail', Fp.filetail);
    end
end

%% Heatraster sorted by mod score
if plotHeatRaster
    figname = 'perCondArea_modSort';
    % need to fix this..
    for a = 1:length(modF)
        animal = modF(a).animal{3};
        Pp=load_plotting_params({'defaults', 'HeatRastSUmod'}, 'pausefigs', pausefigs, ...
            'savefigs', savefigs);
        time = modF(a).output{1}(1).time';
        s = knnsearch(time, -Pp.win(1));
        e = knnsearch(time, Pp.win(2));
        time = time(s:e);
        nrows = length(modF(a).dmatIdx);
        ncols = length(Fp.areas);
        for iv = 1:length(modF(a).dmatIdx) % per condition mean
            for ar = 1:length(Fp.areas) % per area/subarea
                areaIdx = strcmp({modF(a).output{1}.area}', Fp.areas{ar}{1});
                subareaIdx = ~cellfun(@isempty, strfind({modF(a).output{1}.subarea}', Fp.areas{ar}{2}));
                iareaIdx = find(all([areaIdx subareaIdx],2));
                sf = subaxis(nrows,ncols,ar*iv,Pp.posparams{:});
                sf.Tag = strjoin(Fp.areas{:,ar});
                areaMeanStack = cell2mat(arrayfun(@(x) modF(a).output{1}(x).evMean{iv}(s:e), iareaIdx, 'un',0));
                modStack = cell2mat(arrayfun(@(x) modF(a).output{1}(x).mPctChange{iv}, iareaIdx, 'un',0));
                [~, modsortIdx] = sort(modStack, 'descend');
                imagesc(time, 1:length(iareaIdx), areaMeanStack(modsortIdx,:));
                caxis([0 300])
                line([0 0], ylim, 'Color', 'k')
                title(sprintf('%s %s', Fp.areas{ar}{1}, Fp.areas{ar}{2}))
                if ar == 1
                    ylabel('cellNum')
                    xlabel('time (s) from SWR')
                end
            end
        end
        %%
        stit = sprintf('%s %s %s', Fp.Label, animal, Fp.env);
        setSuperAxTitle(stit);
        if pausefigs
            pause;
        end
        if savefigs
            save_figure([pconf.andef{4} '/' Fp.Label '/' animal], stit);
        end
    end
end

%% =============================================================================
%% =============================================================================
%% SU Lick Response
%%
if runLickPhaseSpiking
    Fp = [];
    Fp.Label = 'lickSpikeXC';
    Fp.animals = animals; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
    Fp.filtfunction = 'dfa_lickXCorrSpikes';
    Fp.params = {Fp.Label, 'wtrackdays', 'valid_ntrodes', Fp.filtfunction, 'nonMU_cells'};
    Fp = load_filter_params(Fp);
    Fp.areas = {{'ca1', 'd'}, {'mec', 'supf'}, {'mec', 'deep'}};
    Fp.filetail = sprintf('_%s_%s', Fp.eventType, Fp.env);
    %% FF
    if create_filter
        F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes',...
            Fp.tetfilter, 'excludetime', Fp.timefilter,'iterator',Fp.iterator, 'cells',...
            Fp.cellfilter);
        F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    end
    if run_ff
        F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F), 'un', 1);
        F = runfilter(F);
        save_data(F, [pconf.andef{2} Fp.Label], Fp.filtfunction, ...
            'filetail', Fp.filetail);
    end
    if load_ffdata
        F = load_data([pconf.andef{2} Fp.Label], Fp.filtfunction, Fp.animals, ...
            'filetail', Fp.filetail);
    end
end

%% plot per SU ILPC With that cell's area's heatraster.
if plotClustILPC
    figname = 'suILPC';
    Pp=load_plotting_params({'defaults', figname}, 'pausefigs', pausefigs, ...
        'savefigs', savefigs, 'init_fig', 0);
    binphase = linspace(0, 2*pi, Pp.numBins);
    binphC = binphase(1:end-1) + diff(1:2)/2;
    x = [binphC binphC+max(binphC)];
    for a = 1:length(F)
        animal = F(a).animal{3};
        for ar = 1:length(Fp.areas)
            % find cells in this area
            areaIdx = strcmp({F(a).output{1}.area}', Fp.areas{ar}{1});
            subareaIdx = ~cellfun(@isempty, strfind({F(a).output{1}.subarea}', Fp.areas{ar}{2}));
            iareaIdx = find(all([areaIdx subareaIdx],2));
            
            % make heatraster from all the clusters in this area
            spikeLickPhaseHist = cell2mat(arrayfun(@(x) ...
                histcounts(F(a).output{1}(x).spikeLickPhase, binphase, 'Normalization', ...
                'probability'), iareaIdx, 'un', 0));
            spikeLickPhaseHistShift = circshift(spikeLickPhaseHist, round(size(spikeLickPhaseHist,2)/2), 2);
            [~, imx] = max(spikeLickPhaseHistShift,[],2);
            %             [~, modsortIdx] = sort(max(spikeLickPhase,[],2), 1, 'descend');
            [~, modsortIdx] = sort(imx, 1, 'descend');
            spikeLickPhaseSort = spikeLickPhaseHistShift(modsortIdx,:);
            spikeLickPhaseZ = zscore(smoothdata(spikeLickPhaseSort, 2, 'loess', 6), [],2);
            
            for ic = 1:length(iareaIdx) % for each cluster in this area
                c = iareaIdx(ic);
                init_plot(pausefigs, savefigs, 'position', Pp.position);
                day = F(a).output{1}(c).index(1);
                eps = F(a).output{1}(c).index(4:end);
                nt = F(a).output{1}(c).index(2);
                clust = F(a).output{1}(c).index(3);
                fprintf('%s %d %d %d\n', animal, day, nt, clust);
                
                %% smooth average ILPC
                sf = subaxis(3,1,1, Pp.posparams{:});
                sf.Tag = 'ilpc';
                h = spikeLickPhaseHist(ic,:);
                bar(x, [h h], 'facecolor', [.5 .5 .5], 'facealpha', .5, 'edgealpha', 0)
                hold on;
                hs = smoothdata([h h], 2,'loess', Pp.numBins/2);
                area(x, hs, 'facecolor', 'k', 'linewidth', 2, 'facealpha', .1)
                hold off;
                axis tight
                xticks([])
                ylabel('spike probability');
                
                %% all su in area heatraster
                sf = subaxis(3,1,2:3, Pp.posparams{:});
                sf.Tag = 'heatrast';
                imagesc(x, 1:length(iareaIdx), repmat(spikeLickPhaseZ,1,2));
                caxis(sf,'auto')
                colormap('jet');
                line([0 0], ylim, 'Color', 'k')
                ylabel(sprintf('%s SU num', strjoin(Fp.areas{:,ar})))
                %                 title(sprintf('%s %s', Fp.areas{ar}{1}, Fp.areas{ar}{2}))
                xlabel('lick phase rad')
                %%
                %             allAxesInFigure = findall(gcf,'type','axes');
                %             linkaxes(allAxesInFigure, 'x');
                stit = sprintf('%s %s %d %d %d %s %s', figname, animal, day, nt, clust, Fp.env, strjoin(Fp.areas{:,ar}));
                setSuperAxTitle(stit);
                if pausefigs
                    pause;
                end
                if savefigs
                    save_figure([pconf.andef{4} '/' figname '/' animal], stit);
                end
            end
        end
    end
end
%% plot just the heatraster.. this was a predecessor
if plotILPTHeatRaster
    figname = 'spikeILPC';
    for a = 1:length(F)
        animal = F(a).animal{3};
        Pp=load_plotting_params({'defaults', figname}, 'pausefigs', pausefigs, ...
            'savefigs', savefigs);
        time = F(a).output{1}(1).time';
        s = knnsearch(time, -Pp.win(1));
        e = knnsearch(time, Pp.win(2));
        time = time(s:e);
        dSUIdx = F(a).output{1}.index;
        ncols = length(Fp.areas);
        for ar = 1:length(Fp.areas) % per area/subarea
            
            sf = subaxis(1,ncols,ar,Pp.posparams{:});
            sf.Tag = strjoin(Fp.areas{:,ar});
            %%
            areaIdx = strcmp({F(a).output{1}.area}', Fp.areas{ar}{1});
            subareaIdx = ~cellfun(@isempty, strfind({F(a).output{1}.subarea}', Fp.areas{ar}{2}));
            iareaIdx = find(all([areaIdx subareaIdx],2));
            spikeLickPhase = cell2mat(arrayfun(@(x) ...
                histcounts(F(a).output{1}(x).spikeLickPhase, [pi:pi/36:2*pi]'), ...
                iareaIdx, 'un', 0));
            spikeLickPhase = circshift(spikeLickPhase, round(size(spikeLickPhase,2)/2), 2);
            [maxPct, imx] = max(spikeLickPhase,[],2);
            [~, modsortIdx] = sort(imx, 'ascend');
            spikeLickPhase = zscore(smoothdata(spikeLickPhase, 2, 'loess', 6), [],2);
            
            imagesc([-pi:pi/24:pi]', 1:length(iareaIdx), spikeLickPhase(modsortIdx,:));
            caxis(sf,'auto')
            colormap('jet');
            line([0 0], ylim, 'Color', 'k')
            title(sprintf('%s %s', Fp.areas{ar}{1}, Fp.areas{ar}{2}))
            %%
            if ar == 1
                ylabel('cellNum')
                xlabel('% ILI')
            end
        end
        %%
        stit = sprintf('%s %s %s', figname, animal, Fp.env);
        setSuperAxTitle(stit);
        if pausefigs
            pause;
        end
        if savefigs
            save_figure([pconf.andef{4} '/' Fp.Label '/' animal], stit);
        end
    end
end

%% plot the mrv of the sig phasemod su
if plotCumPolar
    figname = 'SULickILPC_cdfPolar';
    for a = 1:length(F)
        animal = F(a).animal{3};
        Pp = load_plotting_params({'defaults', figname});
        for ar = 1:length(Fp.areas)
            ifig = init_plot(showfigs, 'position', Pp.position);
            % find SU from this Area
            iareaIdx = strcmp({F(a).output{1}.area}', Fp.areas{ar}{1});
            %% cdf SU phasemod
            sf = subaxis(2,1,1,Pp.posparams{:});
            sf.Tag = strjoin(Fp.areas{:,ar});
%             subareaIdx = ~cellfun(@isempty, strfind({F(a).output{1}.subarea}', Fp.areas{ar}{2}));
%             iareaIdx = find(all([iareaIdx subareaIdx],2));

            % shuffle logZ cdf
            phasemodShuf = cell2mat([F(a).output{1}(iareaIdx).phasemodShuf]);
            [h, b] = histcounts(phasemodShuf, Pp.logZBins, 'Normalization','cdf');
            logZticks = b(1:end-1)+diff(b)/2;
            plot(logZticks, h, 'k') 
            hold on;
            % real
            phasemod = [F(a).output{1}(iareaIdx).phasemod];
            [h, b] = histcounts(phasemod, Pp.logZBins, 'Normalization','cdf');
            logZticks = b(1:end-1)+diff(b)/2;
            plot(logZticks, h, 'b') 
            axis tight
            % sig line
            phShSigLen = length(phasemodShuf) * Pp.sigPct;
            phShufSort = sort(phasemodShuf);
            siglin = phShufSort(round(phShSigLen));
            line([siglin siglin], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', 1)

%             ax = gca;
%             ax.YDir = 'reverse'; % to match karalis
            xlabel('log(|MRV|)')
            ylabel('% SU');
            grid on            
            xl = xlim;
            xlim([-4 xl(2)])
            hold off
            
            %% sig phasemod mrv POLAR
            sf = subaxis(2,1,2,Pp.posparams{:});
            sf.Tag = strjoin(Fp.areas{:,ar});
            vecang = [F(a).output{1}(iareaIdx).vecang];
            mrvmag = [F(a).output{1}(iareaIdx).meanMRVmag];
            idxSig = phasemod > siglin;
%             polarhistogram(vecang, 36, 'Normalization', 'pdf', 'facecolor', 'b')
%             hold on
            [thS, thIdx] = sort(vecang(idxSig));
            r = mrvmag(idxSig);
            rS = r(thIdx);
            polarplot(thS, rS, 'ok', 'markersize', 5, 'markerfacecolor', 'k')
            hold off
            Ax = gca; % current axes
%             Ax.ThetaGrid = 'off';
            Ax.RGrid = 'off';
            thetaticks([0 90 180 270])
            Ax.ThetaAxisUnits = 'radians';
            Ax.RTickLabel = [];
            Ax.ThetaAxis.Label.String = 'MRV sigSU';
%             Ax.ThetaTickLabel = [];
%             [x, y] = pol2cart(thS, rS); 
%             fill(x,y, 'r')

            %%
            stit = sprintf('%s %s %s %s', figname, animal, Fp.env, strjoin(Fp.areas{:,ar}));
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure([pconf.andef{4} '/' figname '/' animal], stit);
                F(a).figs.spikePhaseCumPolar{ar} = strsave;
            end
        end
    end
end








