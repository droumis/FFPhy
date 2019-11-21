%{

So far:
    - dfa_eventTrigSpiking (collected psth and eta; singlecellanal)
        - combineEpochs
        - calcSUmod

    - dfa_lickXCorrSpikes (xcorr, excorr, phasemod; singleDayCellAnal)

need to combine the sumod/psth code with the lickXcorrSpikes code..
- or do i?
- ideally dfa_eventTrigSpiking would call calcSUmod within the FF, but for
that i would need to combine the epochs via singleDayCellAnal and therefore
refactor dfa_eventTrigSpiking and calcSUmod a bit..

idk what to do.. i guess i just want some examples of lick triggered
spiking from SU, and then a population summary of the modulation score for
the population of SU.. and also a ETA of all the SU cells.. which i guess
could either be the ETA from dfa_eventTrigSpiking or the xcorr from
dfa_lickXCorrSpikes...
the 'mod' score could be
    - par: Raleigh, nonpar: excorr, phasemod, sumod..
DO:
== remake the lick triggered rasters with dfa_eventTrigSpiking via singleDayCellAnal
-- raster psth, eta per unit..
-- heatraster eta for all the units, per area
-- shuffle vs for [excorr, phasemod, sumod, raleigh] per unit
-- prop sig per area, per mod score (excorr, phasemod, sumod, raleigh).


%}
pconf = paramconfig;
runEventTrigSpiking = 1;
runLickXCorrSpikes = 0;

create_filter = 0;
run_ff = 0;
load_ffdata = 0;

% via dfa_eventTrigSpiking. per cluster raster
plotClustPSTH = 1;

% via dfa_eventTrigSpiking, calcMod. per animal clusts mean psth, stat
calcMod = 0;
loadMod = 0;
plotHeatRaster = 0;
plotStatMod = 0;

% via dfa_lickXCorrSpikes
plotXCPH = 0;
plotILIvIBI = 0;

pausefigs = 0;
savefigs = 1;

%%
if runEventTrigSpiking
    Fp = [];
    Fp.Label = 'lickTrigSUmod';
    Fp.animals = {'JZ1'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
    Fp.filtfunction = 'dfa_eventTrigSpiking';
    Fp.params = {Fp.Label, 'wtrackdays', 'valid_ntrodes', Fp.filtfunction, 'nonMU_cells'};
    Fp = load_filter_params(Fp);
    Fp.areas = [{'ca1 d'}, {'mec supf'}, {'mec deep'}];
    
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
        save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
            'filetail', sprintf('_%s_%s', Fp.env))
    end
    if load_ffdata
        F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
            'filetail', sprintf('_%s_%s', Fp.env));
    end
end
%% plot per SU cluster raster PSTH
if plotClustPSTH
    for a = 1:length(F)
        animal = F(a).animal{3};
        cellInfo = loaddatastruct(F(a).animal{2}, animal, 'cellinfo');
        numClusts = length(F(a).output{1});
        for c = 1:numClusts % per SU
            Pp=load_plotting_params({'defaults', Fp.Label}, 'pausefigs', pausefigs, ...
                'savefigs', savefigs);
            day = F(a).output{1}(c).index(1);
            eps = F(a).output{1}(c).index(4:end);
            nt = F(a).output{1}(c).index(2);
            clust = F(a).output{1}(c).index(3);
            fprintf('%s %d %d %d\n', animal, day, nt, clust);
            
            %% lick triggered raster plot
            sf = subaxis(3,1,[1 2],Pp.posparams{:});
            sf.Tag = 'raster';
            % trim
            time = F(a).output{1}(c).time';
            s = knnsearch(time, -Pp.win(1));
            e = knnsearch(time, Pp.win(2));
            time = time(s:e);
            
            [xx,yy] = find(F(a).output{1}(c).psth(:,s:e)');
            %             spikeClr =
            h = scatter(xx*Fp.bin-Pp.win(1), yy, Pp.spikeSz, '.k', 'markeredgealpha', Pp.spikeAlpha);
            axis tight
            xlim([-Pp.win(1) Pp.win(2)])
            xticks([]);
            ylabel('lickNum');
            line([0 0], ylim, 'linestyle', '--', 'color', [1 .5 .5 .5])
            
            %% smooth average/PSTH
            % spike psth
            sf = subaxis(3,1,3,Pp.posparams{:});
            sf.Tag = 'psth';
            
            time = F(a).output{1}(c).time';
            bintime = [time(1):Pp.bin:time(end)]';
            s = knnsearch(bintime, -Pp.win(1));
            e = knnsearch(bintime, Pp.win(2));
            
            [xx,yy] = find(F(a).output{1}(c).psth');
            h = histc(xx*Fp.bin-Pp.win(1), bintime); %,'Normalization', 'probability', ...
            %                 'edgealpha', 0, 'facecolor', 'k', 'facealpha', 1)+
            hs = smoothdata(h, 1,'loess',10);
            area(bintime(s:e), hs(s:e), 'facecolor', 'k')
            axis tight
            xlabel('time from lick s')
            line([0 0], ylim, 'linestyle', ':', 'color', [1 .5 .5 .5])
            
            %%
            allAxesInFigure = findall(gcf,'type','axes');
            linkaxes(allAxesInFigure, 'x');
            stit = sprintf('%s %s %d %d %d %s', Fp.Label, animal, day, nt, clust, Fp.env);
            setSuperAxTitle(stit);
            
            if pausefigs
                pause;
            end
            if savefigs
                save_figure([pconf.andef{4} '/' Fp.Label '/' animal], stit);
            end
        end
    end
end

%% =============================================================================
%% =============================================================================
%% calc event Trig Modulation
if calcMod
    modF = calcSUmod(F);
end
if loadMod
    modF = load_data('results', 'sumod', Fp.animals);
end
%% Heatraster sorted by mod score
if plotHeatRaster
    figname = 'HeatRastSUmod';
    % need to fix this..
    for a = 1:length(modF)
        animal = modF(a).animal{3};
        if 1 % i need to have added the info to the data structs before this
            cellInfo = loaddatastruct(F(a).animal{2}, animal, 'cellinfo');
            cInfo = cellfetch(cellInfo, '');
        end
        Pp=load_plotting_params({'defaults', figname}, 'pausefigs', pausefigs, ...
            'savefigs', savefigs);
        numClusts = length(mF(a).output{1});
        % get common time
        time = mF(a).output{1}(1).time';
        bintime = [time(1):Pp.bin:time(end)]';
        s = knnsearch(bintime, -Pp.win(1));
        e = knnsearch(bintime, Pp.win(2));
        psth = {};
        suIdx = cell2mat({mF(a).output{1}.index}');
        for ar = 1:length(areas)
            [~, suIdxCInfo] = ismember(suIdx(:,[1 4 2 3]),cInfo.index(:,[1 2 3 4]),'rows');
            numClustsArea
            for c = 1:numClustsArea % per SU, get psth
                day = mF(a).output{1}(c).index(1);
                eps = mF(a).output{1}(c).index(4:end);
                nt = mF(a).output{1}(c).index(2);
                clust = mF(a).output{1}(c).index(3);
                area = cellInfo{day}{eps(1)}{nt}{clust}.area;
                subarea = cellInfo{day}{eps(1)}{nt}{clust}.subarea;
                fprintf('%s %d %d %d %s %s\n', animal, day, nt, clust, area, subarea);
                
                [xx,yy] = find(mF(a).output{1}(c).psth');
                h = histc(xx*Fp.bin-Pp.win(1), bintime);
                hs = smoothdata(h, 1,'loess', 10);
                ar = find(strcmp([area, ' ', subarea], Fp.areas));
                if isempty(ar)
                    continue
                end
                psth{ar} = [psth{ar}; h(s:e)];
            end
            imagesc(time, 1:size(FRrast,1), FRrast(clusterRank,:))
        end
        
        %%
        %         allAxesInFigure = findall(gcf,'type','axes');
        %         linkaxes(allAxesInFigure, 'x');
        stit = sprintf('%s %s %d %d %d %s', figname, animal, day, Fp.env);
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
%% Run dfa_lickXCorrSpikes

%%
if runLickXCorrSpikes
    Fp = [];
    Fp.Label = 'lickSpikeXC';
    Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
    Fp.filtfunction = 'dfa_lickXCorrSpikes';
    Fp.params = {Fp.Label, 'wtrackdays', 'valid_ntrodes', Fp.filtfunction, 'nonMU_cells'};
    Fp = load_filter_params(Fp);
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
        save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
            'filetail', sprintf('_%s_%s', Fp.env))
    end
    if load_ffdata
        F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
            'filetail', sprintf('_%s_%s', Fp.env));
    end
end
%% Plot XC HeatRaster, EXCORR, PolarRaleigh, Phasemod
if plotXCPH
    figname = 'xcorrPhasemodSU';
    for a = 1:length(F)
        animal = F(a).animal{3};
        for c = 1:length(F(a).output{1})
            idata = F(a).output{1}(c);
            day = idata.index(1);
            nt = idata.index(2);
            clust = idata.index(3);
            if isempty(idata.phasemod)
                continue
            end
            if sum(idata.numSpikesPerEp) < 100
                continue
            end
            if sum(idata.numLicksPerEp) < 10
                continue
            end
            if length(idata.spikeLickPhaseShuf) < 10
                continue
            end
            Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
                'savefigs', savefigs);
            
            %% xcorr
            sf = subaxis(2,2,1,Pp.posparams{:});
            sf.Tag = 'xcorr';
            
            % add PPwin trim
            % shuffled xcorr with std ghost trail
            time = idata.xc.time';
            s = knnsearch(time, -Pp.win(1));
            e = knnsearch(time, Pp.win(2));
            time = time(s:e,1);
            xmsh = nanmean(cell2mat(idata.smthxcShuf'))';
            xmsh = xmsh(s:e,1);
            xstdsh = nanstd(cell2mat(idata.smthxcShuf'))'; %/size(cell2mat(out.smthxcShuf'),1);
            xstdsh = xstdsh(s:e,1);
            
            plot(time, xmsh, 'color', [0 0 1 .2], 'linewidth', 1);
            hold on;
            fill([time; flipud(time)],[xmsh-xstdsh;flipud(xmsh+xstdsh)],'b', 'linestyle', ...
                'none', 'facealpha', .2);
            % xcorr norm
            bar(time, idata.normxc(s:e), 'k', 'facealpha', .2, 'edgealpha', 0)
            % xcorr norm smooth
            plot(time, idata.smthxc(s:e), 'k')
            line([0 0], ylim, 'color', 'k', 'linestyle', '--', 'linewidth', .5)
            
            ylabel('xcorr');
            xlabel('time from lick s');
            axis tight
            hold off;
            
            %% excorr
            sf = subaxis(2,2,2,Pp.posparams{:});
            sf.Tag = 'excorr';
            histogram(cell2mat(idata.excorrShuf), 60,'Normalization','probability','edgealpha', 0, 'facecolor', 'k');
            excsort = sort(cell2mat(idata.excorrShuf));
            idxsig = round(Pp.sigpct/100*length(idata.excorrShuf));
            line([excsort(idxsig) excsort(idxsig)], ylim, 'color', [0 0 0 .8], 'linestyle', '--');
            hold on;
            line([idata.excorr idata.excorr], ylim, 'color', 'r');
            excp = 1-sum(idata.excorr>cell2mat(idata.excorrShuf))/length(idata.excorrShuf);
            title(sprintf('excorr %.03f p%.03f', idata.excorr, excp));
            ylabel('probability')
            xlabel('excess corr')
            axis tight
            hold off;
            
            %% phase
            sf = subaxis(2,2,3,Pp.posparams{:});
            sf.Tag = 'phase';
            polarplot([zeros(size(idata.spikeLickPhase,1),1) idata.spikeLickPhase]', ...
                repmat([0 1],size(idata.spikeLickPhase,1),1)', 'color', [0 0 0 .4], 'linewidth', 4);
            hold on
            polarplot([0; idata.vecang], [0; idata.meanMRVmag], 'color', [1 0 .3], 'linewidth', 4)
            grid off
            rticks([])
            thetaticks([])
            title('swr ILI-phase')
            hold off
            axis tight
            [Rp, z] = circ_rtest(idata.spikeLickPhase);
            title(sprintf('Rp%.03f Rz%.03f', Rp, z));
            
            %% phasemod
            sf = subaxis(2,2,4,Pp.posparams{:});
            sf.Tag = 'phasemod';
            histogram(cell2mat(idata.phasemodShuf), 100, 'Normalization', 'probability', 'edgealpha', 0, 'facecolor', 'k');
            hold on;
            mrvsort = sort(cell2mat(idata.phasemodShuf));
            idxsig = round(Pp.sigpct/100*length(idata.phasemodShuf));
            line([mrvsort(idxsig) mrvsort(idxsig)], ylim, 'color', [0 0 0 .8], 'linestyle', '--');
            hold on
            line([idata.phasemod idata.phasemod], ylim, 'color', 'r');
            modp = 1-sum(idata.phasemod>cell2mat(idata.phasemodShuf))/length(idata.phasemodShuf);
            title(sprintf('logMRVmag %.03f p%.03f', idata.phasemod, modp));
            ylabel('probability')
            xlabel('log(Rayleigh Z)')
            axis tight
            hold off
            
           %%
           stit = sprintf('%s %s %d %d %d %s', figname, animal, day, nt, clust, Fp.env);
           setSuperAxTitle(stit);
           
           if pausefigs
               pause;
           end
           if savefigs
               save_figure([pconf.andef{4} '/' Fp.Label '/' animal], stit);
           end
        end
    end
end
%% Plot LM time/pct lick ILI vs Burst IBI
if plotILIvIBI
    figname = 'ILIvIBI';
    for a = 1:length(F)
        animal = F(a).animal{3};
        for c = 1:length(F(a).output{1})
            idata = F(a).output{1}(c);
            day = idata.index(1);
            nt = idata.index(2);
            clust = idata.index(3);
            if isempty(idata.phasemod)
                continue
            end
            if sum(idata.numSpikesPerEp) < 100
                continue
            end
            if sum(idata.numLicksPerEp) < 10
                continue
            end
            if length(idata.spikeLickPhaseShuf) < 10
                continue
            end
            Pp=load_plotting_params({'defaults',figname}, 'pausefigs', pausefigs, ...
                'savefigs', savefigs);
            
            %% lin fit of time 
            sf = subaxis(2,2,1,Pp.posparams{:});
            x = idata.spikeTimeSinceBurstStart;
            y = idata.spikeTimeSinceLick;
            plot(x, y, '.')
            xlabel('time since burst start (s)')
            ylabel('time since lick start (s)')
            
            coefficients = polyfit(x, y, 1);
            xFit = linspace(min(x), max(x), 1000);
            yFit = polyval(coefficients , xFit);
            hold on;
            plot(xFit, yFit, 'r-', 'LineWidth', 2);
            grid on;
            
            %% lin fit of pct
            sf = subaxis(2,2,2,Pp.posparams{:});
            spikePctSinceBurst = idata.spikeTimeSinceBurstStart ./ diff(idata.spikeBurstInterval,[],2);
            x = spikePctSinceBurst;
            y = idata.spikePctSinceLick;
            plot(x, y, '.')
            xlabel('pct IBI (%)')
            ylabel('pct ILI (%)')
            
            coefficients = polyfit(x, y, 1);
            xFit = linspace(min(x), max(x), 1000);
            yFit = polyval(coefficients , xFit);
            hold on;
            plot(xFit, yFit, 'r-', 'LineWidth', 2);
            grid on;
            
           %% lin fit of order
            sf = subaxis(2,2,3,Pp.posparams{:});
            x = idata.spikeBurstLickNum;
            y = idata.spikePctSinceLick;
            plot(x, y, '.')
            xlabel('Burst Lick order')
            ylabel('pct ILI (%)')
            
            coefficients = polyfit(x, y, 1);
            xFit = linspace(min(x), max(x), 1000);
            yFit = polyval(coefficients , xFit);
            hold on;
            plot(xFit, yFit, 'r-', 'LineWidth', 2);
            grid on;
            
            %%
            stit = sprintf('%s %s %d %d %d %s', figname, animal, day, nt, clust, Fp.env);
            setSuperAxTitle(stit);
            
            if pausefigs
                pause;
            end
            if savefigs
                save_figure([pconf.andef{4} '/' Fp.Label '/' animal], stit);
            end
        end
    end
end

