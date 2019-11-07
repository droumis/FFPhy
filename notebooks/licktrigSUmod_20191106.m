
%{
licktrig su
DR 110519

fuck.. merge dfa_licktrigspiking.m into dfa_eventTrigSpiking
clean the calc mod stuff up
%}

% get data
pconf = paramconfig;
create_filter = 0;
run_ff = 0;
load_ffdata = 0;
% stack data
combineEpochs = 0;
loadCombinedEpochs = 0;
% design mat
createDM = 0;
% get event response
calcmod = 0;
loadmodF = 0;
% plot
plotfigs = 0;
plotPSTH = 1;
plotClusterFigs = 0;
plotPopulationFigs = 0;

pausefigs = 1;
savefigs = 0;

%% filter params
Fp.animals = {'D10'}; %, 'D12', 'D13', 'JZ1', 'JZ2', 'JZ4'};
Fp.filtfunction = 'dfa_eventTrigSpiking';
Fp.Label = 'lickTrigSpiking';
Fp.params = {'licks', 'exemplar_wepochs', 'valid_ntrodes', 'nonMU_cells', Fp.filtfunction};
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
        'filetail', ['_' Fp.env '_' Fp.eventType]);
end
if load_ffdata
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, Fp.animals, ...
        'filetail', ['_' Fp.env '_' Fp.eventType]);
end
%% combine data

if combineEpochs
    cF = combine_epochs(F, Fp, Fp.paths);
end
if loadCombinedEpochs
    cF = load_data(Fp.paths.filtOutputDirectory, Fp.Label, Fp.animals,...
        'filetail', '_combEps');
end
%% create per epoch designmat...
if createDM
    % given a list of conditions and the animal, day, epoch, eventTimes, create the design mat
    for a = 1:length(F)
        data(a).day = [];
        data(a).epoch = [];
        data(a).evStart = [];
        day = 0;
        epoch = 0;
        for c = 1:length(F(a).output{1})
            % 1 per unique epoch
            if F(a).output{1}{c}.index(1:2) == [day epoch]
                continue
            end
            day =  F(a).output{1}{c}.index(1);
            epoch = F(a).output{1}{c}.index(2);
%             try
                data(a).animal = F(a).animal{3};
%             catch
%                 data(a).animal = F(a).animal;
%             end
            numEvents = numel(F(a).output{1}{c}.eventTimes);
            data(a).day = [data(a).day repmat(day,numEvents,1)];
            data(a).epoch = [data(a).epoch repmat(epoch,numEvents,1)];
            data(a).evStart = [data(a).evStart F(a).output{1}{c}.eventTimes];
        end
    end
    defaults = {'exemplar_wepochs'};
    dmat = makeExpvarCatDesignMat(data, 'expvars', {'all', 'lickbouts', 'nolickbouts'}, ...
        'defaults', defaults);
end

%% ---------------- calc mod---------------------------------------
if calcmod
    f = struct;
    f.params = {'wtrackdays', 'lickbouts', 'ca1SU', 'lickTrigSpikingMod'};
    f = load_filter_params(f);
    for a = 1:length(cF)
%         f.animal = ppF.animal;
%         cF = createfilter('animal', f.animal, 'epochs', f.epochfilter,...
%             'excludetime', f.timefilter, 'cells', f.cellfilter);
        modF = calcSUmod(dmat, cF, 'respwin', f.respwin, 'basewin', f.basewin, 'minNumSwr', ...
            f.minNumSwr, 'nshuffs', f.nshuffs, 'shuffms', f.shuffms, 'filetail', ['_' f.env]);
    end
end
if loadmodF
    modF = load_data([pconf.andef{3} '/sumod'], 'sumod', Fp.animals, 'filetail', ['_' Fp.env]);
end
%% ---------------- plot ---------------------------------------

%% plot lick trig mean heatrasters sorted by mod
if plotfigs
    if plotPSTH
        figname = 'lickTrigSpikingPSTH';
        for a = 1:length(cF)
            animal = cF.animal;
            for ic = 1:length(cF(a).data)
                Pp=load_plotting_params({'defaults',figname});
                %% raster
                sf1 = subaxis(2,1,1,Pp.posparams{:});
                set(gca, 'Tag', 'raster');
                h=zoom; h.Motion = 'horizontal'; h.Enable = 'on'; p=pan; p.Motion='horizontal';
                psth = cF(a).data(ic).psth;
                [xx, yy] = find(psth');
                f = scatter(xx/1000-1.001,yy, Pp.psthSize, 'filled', 'MarkerFaceColor', 'k');
                f.Marker = 'd';
                f.MarkerEdgeAlpha = 0;
                f.MarkerFaceAlpha = 0.3;
                set(gca, 'TickDir','out', 'Xticklabel',[]);
                axis tight;
                ylabel('licknum','FontSize',12,'FontWeight','bold', 'FontName','Arial')
                line([0 0],ylim, 'color','red', 'linewidth',2, 'Color', [1 0 0 .5]);
                
              %% psth
                sf2 = subaxis(2,1,2);
                set(gca, 'Tag', 'psth');
                h = histogram(xx/1000-1.001,cF(a).data(ic).time(1:20:end), 'facecolor', 'k', 'Normalization', ...
                    'pdf');
                h.EdgeColor = 'none';
                axis tight
                line([0 0],ylim, 'color','red', 'linewidth',2, 'Color', [1 0 0 .5]);
                ylabel('pdf','FontSize',12,'FontWeight','bold', 'FontName','Arial')
                set(gca, 'YGrid', 'off', 'XGrid', 'on','TickDir','out', 'TickLength', [0.001 0]);
                xticklabels([])
                hold off;
                %%
                
                %% ----- link x axis -----
                allAxesInFigure = findall(gcf,'type','axes'); 
                linkaxes(allAxesInFigure, 'x');
                setSuperAxTitle([figname ' ' animal]);
                if pausefigs
                    pause;
                end
                if savefigs
                    save_figure([pconf.andef{4} '/' Fp.Label '/'], [figname animal]);
                end
            end
        end
    end
end

%% plot individual cells lick trig mean +- std/ 

%% plot lick-spike xcorr, excess corr v shuf

%% plot ili-spike phase preference, logMagMRV v shuf

%% per condition per area pct of mod cells



