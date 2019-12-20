
%{
- *** valid_ntrodes filter is not currently being respected for spiking????

- need to finish XP mod SU figure
.. what is the actual final result?
.. where is the result saved? is there an accessible record of per cell
results?
.. what is the final result of the swr mod su?
-- behavior position fig..
%}

pconf = paramconfig;
eventTrigLFP = 0; % PIPE:forest.bear.cactus.mushroom.beer.leaf == eventSet mean Spect
eventTrigSpiking = 1; % PIPE:barn.rat.beer.wheelbarrow == eventSet SU mod
eventType = 'swr'; %lick swr

create_filter = 0;
run_ff = 0;
load_ffdata = 1;

make_dmat = 0;
load_dmat = 0;

calcSUphaseMod = 0;
loadSUphaseMod = 0;

plotPhaseMod_pClust = 0;
plotPhaseMod_HeatRast_pAn = 0;
plotPhaseMod_HeatRast_AllAn = 0;
plotPhaseMod_stats = 0;

calcSUtimeMod = 0;
loadSUtimeMod = 0;

plotTimeMod_pClust = 1;
plotTimeMod_HeatRast_pAn = 0;
plotTimeMod_HeatRast_AllAn = 0;
plotTimeMod_stats = 0;

savefigs = 1;
pausefigs = 0;
showfigs = 0;
savefigas = {'png', 'eps'};

% data filter params
Fp.animals = {'JZ1'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'JZ2', 'JZ4'};

if eventTrigLFP
    Fp.filtfunction = 'dfa_eventTrigLFP'; % Bellicose Bear
    if strcmp(eventType, 'lick')
        expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
        Fp.Label = 'wtrackLickTrigLFP';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
        'excludeAfterLastWell', 'referenced', '4-350Hz', 'lickboutlicks', Fp.Label, Fp.filtfunction};
    elseif strcmp(eventType, 'swr')
        expvars = {'all', 'lickbouts', 'nolickbouts'};
        Fp.Label = 'wtrackSWRTrigLFP';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
        'excludeAfterLastWell', 'referenced', '4-350Hz',  'ripples', ...
        Fp.Label, Fp.filtfunction};
    end
    
    
elseif eventTrigSpiking
    Fp.filtfunction = 'dfa_eventTrigSpiking'; % Redolent Rat
    if strcmp(eventType, 'lick')
        expvars = {'all', 'wetLickBursts', 'dryLickBursts'};
        Fp.Label = 'wtrackLickTrigSpiking';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
            'excludeAfterLastWell', 'nonMU_cells', Fp.Label, Fp.filtfunction};
    elseif strcmp(eventType, 'swr') %'excludeNoise', 
        expvars = {'all', 'lickbouts', 'nolickbouts'};
        Fp.Label = 'wtrackSWRTrigSpiking';
        Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
            'excludeAfterLastWell', 'nonMU_cells', 'ripples', ...
            Fp.Label, Fp.filtfunction}; % 'excludeNoise',
    end
end


% Fp.filtfunction = 'dfa_eventTrigSpiking';
% Fp.Label = 'wtrackLickTrigSpiking';
% Fp.params = {'wtrackdays', 'valid_ntrodes', 'excludePriorFirstWell', ...
%     'excludeAfterLastWell', 'nonMU_cells', Fp.Label, Fp.filtfunction};

%% FF
Fp = load_filter_params(Fp);
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes',...
        Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator, 'cells',...
        Fp.cellFilter);
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

%% make design mat to slice trials
if make_dmat
    dmat = makeExpvarCatDesignMat(F, Fp.expvars, 'eventType', Fp.eventType);
end

if load_dmat
    outdir = 'expvarCat';
    outpath = [pconf.andef{2},outdir,'/'];
    dmat = load_data(outpath, [outdir,'_',Fp.env,'_',Fp.eventType], Fp.animals);
end

%% Spike Phase mod
if calcSUphaseMod % saw
    pmodF = calcPhaseMod(F, dmat);
    save_data(pmodF, 'results', [Fp.Label '_phasemod']);
end
if loadSUphaseMod
    pmodF = load_data('results', [Fp.Label '_phasemod'], Fp.animals);
end

%% calc su modulation
if calcSUtimeMod % wheelbarrow
    modF = calcSUmod(F, dmat);
    save_data(modF, 'results', [Fp.Label '_timemod']);
end
if loadSUtimeMod
    modF = load_data('results', [Fp.Label '_timemod'], Fp.animals);
end
    

%% PLOT=====================================================================

%% Phase Mod PER CLUST ( + per dmc)
if plotPhaseMod_pClust
    % plot the phase mod as a circular spiking probability of IEI elapsed
    figname = 'PhaseModperClust';
    Pp=load_plotting_params({'defaults', figname});
    % create cell array (an x dmc) of struct array perC
    f = reshape([pmodF.output]', length(pmodF(1).output), [])';
    allpmod = arrayfun(@(x) cat(1,f{:,x}), 1:size(f,2), 'un', 0);
    binphaseEdges = linspace(0, 2*pi, Pp.numBins);
    binphCenters = binphaseEdges(1:end-1) + diff(1:2)/2;
    for c = 1:length(allpmod{1})
        ifig = init_plot(showfigs, Pp.position); % fig per clust
        for idmc = 1:size(allpmod,2)
            sf = subaxis(1,size(allpmod,2),idmc,Pp.posparams{:}); % subfig per dmat condition
            cD = allpmod{idmc}(c);
            % get polar prob hist
            spikeIEIphase = cD.spikeIEIPhase;
            spikeIEIPhaseHistProb = histcounts(spikeIEIphase, binphaseEdges, ...
                'Normalization', 'probability');
            if strcmp(pmodF(1).outputIdx{idmc}, 'all')               
                rlim = ceil(max(spikeIEIPhaseHistProb)*100)/100;
                clarea = cD.area;
                clsubarea = cD.subarea;
            end
            hold on;
            % draw polar axes on cartesian coordinatess
            rlim = ceil(max(spikeIEIPhaseHistProb)*100)/100;
            rrv = [0:90:360]*pi/180; % Calculate Radials
            [rd00x, rd00y] = pol2cart(rrv, zeros(size(rrv)));
            [rd15x, rd15y] = pol2cart(rrv, rlim*ones(size(rrv)));
            rdx = [rd00x; rd15x];
            rdy = [rd00y; rd15y];
%             cfull = linspace(0, 2*pi, 1000); % Full Circle Circumference
%             ccf = @(rd,th) [rd*cos(th); rd*sin(th)]; % Function
%             cc1 = ccf(rlim,cfull);
%             cc2 = ccf(rlim*(2/3),cfull);
%             cc3 = ccf(rlim*(1/3),cfull);
            cr = .1;
%             plot(cc1(1,:), cc1(2,:), 'color', [cr cr cr], 'linewidth', 1) % Plot Full Circumference
%             plot(cc2(1,:), cc2(2,:), 'color', [cr cr cr], 'linewidth', 1)
%             plot(cc3(1,:), cc3(2,:), 'color', [cr cr cr], 'linewidth', 1)
            plot(rdx, rdy, 'color', [cr cr cr], 'linewidth', 1) % Plot Radials
            % convert pol2cart
            [x1, y1] = pol2cart(binphCenters, smoothdata(spikeIEIPhaseHistProb));
            % plot smooth polar probability
            h1 = fill(x1, y1, 'k', 'facecolor', [.4 .4 .4]);
            
            set(h1,'facealpha',.1,'edgecolor','k', 'linewidth', 2);
            set(gca, 'Box','off', 'XColor','none', 'YColor', 'none')
            xlim([-(max(abs(xlim)+.01)) max(abs(xlim))+.01])
            axis equal
            xticks([])
            yticks([])
            pbaspect(sf, [1 1 1])
            title(sprintf('dm %s', pmodF(1).outputIdx{idmc}))
            % add result info as text box in lower left of subfig
            pos = get(sf, 'position');
            dim = pos.*[1 1 0.5 0.5];
            annotation('textbox', dim, 'String', {['p = ' num2str(cD.pval)],...
                ['% rank = ', num2str(cD.modPctRank)],['phmod = ', num2str(cD.phasemod)],...
                ['nSpk = ', num2str(length(cD.spikeIEIPhase))],...
                ['nIEI = ', num2str(length(cD.ILI))]},...
                'FontSize', 14, 'FontName','Arial', 'LineStyle','none', 'vert', 'bottom', 'FitBoxToText','on');
            hold off
        end
        stit = sprintf('%s %d %d %d %s %s %s', cD.animal, cD.index(1), cD.index(2), cD.index(3),...
            clarea, clsubarea, figname);
        setSuperAxTitle(stit);
        if pausefigs
            pause
        end
        if savefigs
            strsave = save_figure([pconf.andef{4} '/' figname '/' cD.animal], stit, 'savefigas', savefigas);
        end
    end
end

%% Phase Mod HeatRast PER ANIMAL ( + per dmc, area)
if plotPhaseMod_HeatRast_pAn

    areaLab = cellfun(@(x) [{x{idmc}.area}' {x{idmc}.subarea}'], {pmodF(a).output}', 'un', 0); % per animal [area subarea]
    areaIdx = all([strcmp(areaLab{:}, Fp.areas{ir}{1}) + strcmp(areaLab{:}, Fp.areas{ir}{2})],2); % per animal, per area
    
    pvar = cell2mat(cellfun(@(x) [x{idmc}(areaIdx).phasemod], {pmodF(a).output}, 'un', 0)); % per dmc, per animal, per area
    pvar(isnan(pvar)) = [];
    
end

%% Phase Mod HeatRast ALL ANIMALS ( + per dmc, per area)

if plotPhaseMod_HeatRast_AllAn
    figname = 'wXPphaseModSU_HeatRast_AllAn';
    Pp=load_plotting_params({'defaults', figname});
    binphaseEdges = linspace(0, 2*pi, Pp.numBins);
    binphCenters = binphaseEdges(1:end-1) + diff(1:2)/2;
%     x = [binphC]; % binphC+max(binphC)];
    
    for idmc = 1:length(pmodF(1).output) % for each dmat condition
        for ir = 1:length(Fp.areas) % for each area
            pvarAll = {};
            spikeLickPhaseHist = {};
            for a = 1:length(pmodF)
                areaLab = cellfun(@(x) [{x{idmc}.area}' {x{idmc}.subarea}'], {pmodF(a).output}', 'un', 0); % per animal [area subarea]        
                areaIdx = all([strcmp(areaLab{:}, Fp.areas{ir}{1}) + strcmp(areaLab{:}, Fp.areas{ir}{2})],2); % per animal, per area
                pvar = cell2mat(cellfun(@(x) [x{idmc}(areaIdx).phasemod], {pmodF(a).output}, 'un', 0)); % All An, per dmc, per area
                pvar(isnan(pv)) = [];
                pvarAll{a} = pvar;
                for c = 1:sum(areaIdx)
                    spikeLickPhaseHist{a}{c} = histcounts(pmodF(a).output{idmc}(areaIdx(c)).spikeLickPhase,...
                        binphaseEdges, 'Normalization', 'probability');
                end
            end
            d = spikeLickPhaseHist{a}{:};
            spikeLickPhaseHistAll = nanmean(d{:});
            
            pvarAll = cell2mat(pvarAll);
            
            [~, imx] = max(spikeLickPhaseHist,[],2);
            %             [~, modsortIdx] = sort(max(spikeLickPhase,[],2), 1, 'descend');
            [~, modsortIdx] = sort(imx, 1, 'descend');
            spikeLickPhaseSort = spikeLickPhaseHist(modsortIdx,:);
            
            % plot all animals
            ifig = init_plot(showfigs, Pp.position);
            
            % make heatraster from all the clusters in this area,
            % condition
            
            
            
            %                 spikeLickPhaseHistShift = circshift(spikeLickPhaseHist, round(size(spikeLickPhaseHist,2)/2), 2);
            
            spikeLickPhaseZ = zscore(smoothdata(spikeLickPhaseSort, 2, 'loess', 6), [],2);

        end
    end
end

%% Phase Mod stats
if plotPhaseMod_stats
    
end


%% time mod per cluster ( + per dmc)
if plotTimeMod_pClust
    if strcmp(eventType, 'swr')
        figname = 'wtrackSWRSU';
    elseif strcmp(eventType, 'lick')
        figname = 'wtrackLickSU';
    end
    Pp = load_plotting_params({'defaults', figname});
    for a = 1:length(F) % per animal
        animal = F(a).animal{3};
        fulltime = F(a).output{1}(1).time';
        sT = knnsearch(fulltime, -Pp.win(1));
        eT = knnsearch(fulltime, Pp.win(2));
        time = fulltime(sT:eT);
        bintime = [fulltime(1):Pp.bin:fulltime(end)]';
        sB = knnsearch(bintime, -Pp.win(1));
        eB = knnsearch(bintime, Pp.win(2));
        for c = 1:length(F(a).output{1}) % per su
            ifig = init_plot(showfigs, Pp.position);
            day = F(a).output{1}(c).index(1);
            eps = F(a).output{1}(c).index(4:end);
            nt = F(a).output{1}(c).index(2);
            clust = F(a).output{1}(c).index(3);
            carea = F(a).output{1}(c).area;
            csubarea = F(a).output{1}(c).subarea;
            fprintf('%s %d %d %d\n', animal, day, nt, clust);
            
            %% spike raster
            sf = subaxis(3, 1, [1 2], Pp.posparams{:});
            sf.Tag = 'raster';
            [xx,yy] = find(F(a).output{1}(c).psth(:,sT:eT)');
            spikeTwin = (xx*1e-3)-Pp.win(1);
            h = scatter(spikeTwin, yy, Pp.spikeSz, '.k', 'markeredgealpha', ...
                Pp.spikeAlpha);
            axis tight
            xlim([-Pp.win(1) Pp.win(2)])
            xticks([]);
            ylabel('Event #');
            line([0 0], ylim, 'linestyle', '-', 'color', [.5 .5 1 .5], 'linewidth', 2)
            
            %% PSTH
            sf = subaxis(3,1,3,Pp.posparams{:});
            sf.Tag = 'psth';
            [xx,yy] = find(F(a).output{1}(c).psth');
            spikeTwin = (xx*1e-3)-abs(fulltime(1));
            h = histc(spikeTwin, bintime);
            hs = smoothdata(h, 1,'loess', 20);
            area(bintime(sB:eB), hs(sB:eB), 'facecolor', 'k')
            axis tight
            xlabel('time s')
            ylabel('count')
            line([0 0], ylim, 'linestyle', '-', 'color', [.5 .5 1 .5], 'linewidth', 2)
            
            %%
            allAxesInFigure = findall(gcf,'type','axes');
            linkaxes(allAxesInFigure, 'x');
            stit = sprintf('%s %s %d %d %d %s %s', figname, animal, day, ...
                nt, clust, carea, csubarea);
            setSuperAxTitle(stit);
            if pausefigs
                pause
            end
            if savefigs
                strsave = save_figure([pconf.andef{4} '/' figname '/' animal], stit, 'savefigas', savefigas);
            end
        end
    end
end

%% Time Mod PER ANIMAL ( + per dmc, area)
if plotTimeMod_HeatRast_pAn
    if strcmp(eventType, 'swr')
        figname = 'wtrackSWRHeatrast';
    elseif strcmp(eventType, 'lick')
        figname = 'wtrackLickHeatRast';
    end
    Pp = load_plotting_params({'defaults', figname}); % load params
    for a = 1:length(modF) % per animal
        animal = modF(a).animal{3};
        for idmc = 1:length(modF(1).output) % for each dmat condition
            for ir = 1:length(Fp.areas) % per area
                areaLab = cellfun(@(x) [{x{idmc}.area}' {x{idmc}.subarea}'], {modF(a).output}', 'un', 0); % per animal [area subarea]
                areaIdx = all([strcmp(areaLab{:}, Fp.areas{ir}{1}) + strcmp(areaLab{:}, Fp.areas{ir}{2})],2); % per animal, per area
                
                
                pvar = cell2mat(cellfun(@(x) [x{idmc}(areaIdx).mPctChange], {modF(a).output}, 'un', 0)); % per dmc, per animal, per area
                pvar(isnan(pvar)) = [];
            end
        end
    end
end

%% Time Mod HeatRast per animal ( + per dmc, per area)
if plotEventTimeModHeatRast_perAnimal
    if strcmp(Fp.eventType, 'swr')
        figname = 'wtrackSWRHeatrast';
    elseif strcmp(Fp.eventType, 'lick')
        figname = 'wtrackLickHeatRast';
    end
    Pp = load_plotting_params({'defaults', figname}); % load params
    for a = 1:length(modF) % per animal
        animal = modF(a).animal{3};
        for ar = 1:length(Fp.areas) % per area
            numESet = length(modF(a).dmatIdx);
            % find cells in this area
            areaIdx = strcmp({modF(a).output{1}.area}', Fp.areas{ar}{1});
            subareaIdx = ~cellfun(@isempty, strfind({modF(a).output{1}.subarea}', ...
                Fp.areas{ar}{2}));
            iareaIdx = find(all([areaIdx subareaIdx],2));
            iareaIdx = iareaIdx(~arrayfun(@(x) isempty(modF(a).output{1}(x).mPctChange), ...
                iareaIdx,'un',1));
            if isempty(iareaIdx)
                continue
            end
            ifig = init_plot(showfigs, Pp.position); % init fig per area
            for iv = 1:numESet % per eventSet in dmat
                setID = modF(a).dmatIdx{iv};
                sf = subaxis(1, numESet, iv, Pp.posparams{:});
                sf.Tag = 'heatraster';
                mPctChange = [];
                dumpy = [];
                gud = [];
                for i = 1:length(iareaIdx)
                    try
                        if ~isempty(modF(a).output{1}(iareaIdx(i)).mPctChange{iv})
                            
                            mPctChange = [mPctChange; modF(a).output{1}(iareaIdx(i)).mPctChange{iv}];
                            gud = [gud; i];
                            
                        else
                            dumpy = [dumpy; i];
                        end
                    catch
                        dumpy = [dumpy; i];
                    end
                end
                iareaIdx = iareaIdx(gud);
                % make firing rate heatraster from all the clusters in this area
                iFRHRz = cell2mat(arrayfun(@(x) modF(a).output{1}(x).evMeanZ{iv}, ...
                    iareaIdx,'un',0));
                iFRHRzsm = smoothdata(iFRHRz,2,'loess', 10);
                
                %                     pctChange = cell2mat(arrayfun(@(x) modF(a).output{1}(x).mPctChange{iv}, ...
                %                         iareaIdx,'un',0));
                [~, srtIdx] = sort(mPctChange, 1, 'descend');
                iFRHRsmzSorted =  iFRHRzsm(srtIdx,:);
                time = modF(a).output{1}(iareaIdx(1)).time';
                s = knnsearch(time, -Pp.win(1));
                e = knnsearch(time, Pp.win(2));
                time = time(s:e);
                imagesc(time, 1:length(iareaIdx), iFRHRsmzSorted(:,s:e));
                
                colormap(viridis);
                h = colorbar;
                caxis(sf,[-3 3])
                ylabel(h, 'zscore firing rate')
                line([0 0], ylim, 'Color', 'k', 'linestyle', '--')
                title(modF(a).dmatIdx{iv})
                xlabel('time s')
                ylabel('Unit #')
            end
            stit = sprintf('%s %s %s %s', figname, animal, strjoin(Fp.areas{:,ar}));
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

%% Time Mod stats
if plotTimeMod_stats
    
end





