
run_ff = 0;
savedata = run_ff;

load_lfp = 0;
stack_lfp = 0;
calculateAnalyticSignal = 1;
makeStateFilters = 1;
computePWR = 1;
runPermutationTest = 0;
loadPWR = 0;

plotfigs = 0;
savefigs = 1;
pausefigs = 0;

use_filters = {'firstwell', 'noise'};
animals = {'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
add_params = {'wtrack', 'wavelets4-300Hz'};
me = animaldef('Demetris');

Fp.animals = animals;
Fp.add_params = add_params;
Fp.filtfunction = 'dfa_riptriglfp';
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
filetail = '';

%% run filter/func
if run_ff == 1
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, 'eegtetrodes', ...
        Fp.tetfilter, 'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    F = runfilter(F);
    for d = 1:length(F); F(d).datafilter_params = Fp; end
end
%% save data
if savedata == 1
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s%s', Fp.epochEnvironment, filetail))
end
%% load lfp
if load_lfp
    F = load_filter_output(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s%s', Fp.epochEnvironment, filetail));
end
%% stack LFP
if stack_lfp
    lfpstack = stack_riptriglfp(F);
    clear F
end
%% Calculate and Save Analytic Signal
if calculateAnalyticSignal
        computeAnalyticSignal(lfpstack(ian).data{eegidx}, lfpstack(ian).ntrodes, ...
            'saveAnalyticSignal', 1, 'waveSet', Fp.waveSet);
end
%% make statesets as sets of indices into AS ripple dim for specified states
if makeStateFilters
    ripstate = getStateFilters(lfpstack);
end
%% Compute Power from analytic signal for each ntrode, condition
if computePWR
    pwr = getPower(ripstate, Fp, 'baseline', [1 1.5]*1500);
end
%% load power
if loadPWR
    pwr = {};
    for ian = 1:length(Fp.animals)
        savestr = sprintf('%s/power/power_%s_%s_waveSet-%s.mat', me{2}, ...
            Fp.animals{ian}, 'wtrack', Fp.waveSet);
        tmp = load(savestr);
        pwr = [tmp.pwrout; pwr];
    end
end
%% permutation testing two-condition diff
% run this inside the getPower func instead when each AS is loaded..
% if runPermutationTest
%     aIdx = find(ripstate(ian).statesets(:,2)); % rewarded (2)
%     bIdx = find(ripstate(ian).statesets(:,3)); % unrewarded (3)
%     pwr = powerpermtest(pwr, aIdx, bIdx);
% end
%% plotting
% plot each ntrode per condition
% for the difference plots, add the z mask contours
if plotfigs
    Pp = load_plotting_params({'defaults', 'power'});
    for ian = 1:numel(Fp.animals)
        wp = getWaveParams('4-300HzJustfreqs', []);
        animal = Fp.animals{ian};
        anidx = find(strcmp({pwr.animal}, animal));
        for iset = 1:length(pwr(ian).ripstate(1).statesetsfields)
            if savefigs && ~pausefigs
                close all
                ifig =figure('Visible','off','units','normalized','position',Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            for nt = 1:length(pwr(ian).mediandbpower)
                sf = subaxis(5,6,nt, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                    'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, ...
                    'MarginBottom', Pp.MgBm);
                
                idata2plot = squeeze(pwr(ian).mediandbpower{nt}{iset})';
                idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin);
                time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                contourf(sf, time, wp.frex, idata2plot, Pp.contourRes,'linecolor','none');
                set(gca,'ydir','normal','yscale','log'); %,'xlim',[-Pp.pwin(1) Pp.pwin(2)], 'ylim',[wp.frex(1) wp.frex(end)])0                coloraxis = [-5 50]; %'auto';%
                coloraxis = [-1 1]; %[-5 50];
                caxis(coloraxis)
                colormap(Pp.usecolormap)

                hold on;
                title(sprintf('nt%d',nt), 'FontSize',Pp.FontS,'FontWeight',Pp.FontW, 'FontName', ...
                    Pp.FontNm)
                if nt == 1
                    xlabel('time s', 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                    ylabel('freq Hz','FontSize',Pp.FontS, 'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                    ytickskip = 2:4:wp.numfrex;
                    set(gca,'ytick', round(wp.frex(ytickskip)),'yticklabel',round(wp.frex(ytickskip)))
                else
                    set(gca, 'xtick', [])
                    set(gca, 'ytick', [])
                end
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s dbPower %s %s %s', animal, pwr(ian).ripstate(ian).statesetsfields{iset}, add_params{1}, add_params{2});
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 12);
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                me = animaldef('Demetris');
                save_figure(sprintf('%s/wavepower/',me{4}), 'wavepower', sprtit)
                close all
            end
            close all;
        end
    end
end
