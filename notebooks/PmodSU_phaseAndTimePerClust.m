
%{
Replicate phase mod.. include raw spiking in polar plots
add time mod to the same plot as phase mod so it's easier to survey
indivudal units across views
%}

% query
create_filter = 1;
run_ff = 0;
load_ffdata = 0;

% Analysis
make_dmat = 0;
load_dmat = 0;

calcSUphaseMod = 0;
loadSUphaseMod = 0;
calcSUtimeMod = 0;
loadSUtimeMod = 0;

% Plot
plotSUmod_pClust = 0;
plotPhaseMod_pClust = 0;
plotTimeMod_pClust = 0;

savefigs = 0;
pausefigs = 0;
showfigs = 0;
savefigas = {'png','pdf'};

%% Define Filter Params
pconf = paramconfig('Demetris'); % globals per user
Fp.animals = {'JZ4'};
eventType = 'lick';
Fp.filtfunction = 'dfa_eventTrigSpiking';
Fp.Label = 'wtrackLickTrigSpiking'; % used for filename/plots
Fp.params = {'exemplar_wepochs', ...
    'valid_ntrodes', ...
    'excludePriorFirstWell', ...
    'excludeAfterLastWell', ...
    'nonMU_cells', ...
    Fp.Label, ...
    Fp.filtfunction};
Fp = load_filter_params(Fp);

%% Create Filter and Run FiltFunc
if create_filter
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, ...
        'eegtetrodes', Fp.tetfilter, 'excludetime', Fp.timefilter, ...
        'iterator', Fp.iterator, 'cells', Fp.cellFilter);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
end
if run_ff
    F = arrayfun(@(x) setfield(F(x),'datafilter_params',Fp),1:length(F),'un',1);
    F = runfilter(F);
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        'filetail', ['_' Fp.Label]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', ['_' Fp.Label]);
end

%% Create Design Mat
if make_dmat
    dmat = makeExpvarCatDesignMat(F, Fp.expvars, 'eventType', Fp.eventType);
end

if load_dmat
    outdir = 'expvarCat';
    outpath = [pconf.andef{2},outdir,'/'];
    dmat = load_data(outpath, [outdir,'_',Fp.env,'_',Fp.eventType], Fp.animals);
end

%% Event Phase Mod Spikes
if calcSUphaseMod % saw
    pmodF = calcPhaseMod(F, dmat);
    save_data(pmodF, 'results', [Fp.Label '_phasemod']);
end 
if loadSUphaseMod
    pmodF = load_data('results', [Fp.Label '_phasemod'], Fp.animals);
end

%% Event Time Mod Spikes
if calcSUtimeMod % wheelbarrow
    modF = calcSUmod(F, dmat);
    save_data(modF, 'results', [Fp.Label '_timemod']);
end
if loadSUtimeMod
    modF = load_data('results', [Fp.Label '_timemod'], Fp.animals);
end

 %% Plot P and SWR mod per unit
% if plotSUmod_pClust
%     % SWR time mod
%     
%     
%     % P event time mod
%     
%     
%     % polar raster
%     spikeTimes =
%     boutIntv =
%     spikesinbouts = spikeTimes(logical(isExcluded(spikeTimes, boutIntv(1,:))));
%     [N,~,spbin] = histcounts(spikesinbouts, lickTimes);
%     relspiketime = spikesinbouts - lickTimes(spbin);
%     polar([zeros(size(sprads,1),1) sprads]',repmat([0 1],size(sprads,1),1)', 'k')
%     
%     % polar hist
%     
% 
% end
%% Plot Phase Mod PER CLUST ( + per dmc)
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

%% Plot time mod per cluster ( + per dmc)
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
