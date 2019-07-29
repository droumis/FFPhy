load_stack = 0;

get_ripstate = 0;
load_ripstate = 0;
get_designmat = 0;
load_designmat = 1;
load_raw_pwr = 1;

run_regression = 1;
load_regression = 0;

plot_regression = 1;
pausefigs = 0;
savefigs = 1;

run_TFVarcorr = 0;
plot_corrheatmap = 0;
run_multiple_regress = 0;

Fp.animals = {'JZ3'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp.add_params = {'wtrackdays', 'excludeNoise','excludePriorFirstWell', '<4cm/s', ...
    'wavelets4-300Hz'};

Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
Fp.uselfptype = 'eeg';
Fp.useripstates = {'onlywdays','rewarded', 'unrewarded', 'inbound' , 'outbound'};
pconf = paramconfig;
me = animaldef('Demetris');

%% load lfpstack
if load_stack
    lfpstack = load_data(Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack', Fp.animals);
end

%% make, get ripstates.. 

if get_ripstate
    ripstate = getStateFilters(lfpstack);
    save_data(ripstate, [pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack');
end
if load_ripstate
    ripstate = load_data([pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack', Fp.animals);
end
%% make design matrix
if get_designmat
    designmat = getDesignMatrix(lfpstack);
    save_data(designmat, [pconf.andef{2}, 'designmat/'], 'designmat_wtrack');
end
if load_designmat
    designmat = load_data([pconf.andef{2}, 'designmat/'], 'designmat_wtrack', Fp.animals);
end

%% load raw power
if load_raw_pwr
    wp = getWaveParams(Fp.waveSet);
    rawpwr = load_data(sprintf('%s/analyticSignal/', me{2}), ...
        sprintf('AS_waveSet-%s_%s_power', wp.waveSet, Fp.uselfptype), Fp.animals);
end
%% Multi regression
if run_regression
    PV = runDesignDataRegression(designmat, rawpwr);
end
if load_regression
    PV = load_data([pconf.andef{2}, 'powerVarCorr/'], 'powerVarCorr_wtrack', Fp.animals);
end

%% plot regression
if plot_regression
    Pp = load_plotting_params({'defaults', 'powerVarCorr'});
    wp = getWaveParams('4-300Hz');
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        for iv = 1:length(PV(ian).expvars)
            if savefigs && ~pausefigs
                close all
                ifig =figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            numrows = 2;
            numcols = ceil(length(ntrodes) / 2);            
            for nti = 1:length(ntrodes)
                sf = subaxis(numrows,numcols,nti, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                    'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                    Pp.MgTp, 'MarginBottom', Pp.MgBm);
                nt = ntrodes(nti);
                area = ntinfo{1}{1}{nt}.area;
                subarea = ntinfo{1}{1}{nt}.subarea;
                if isnumeric(subarea)
                    subarea = num2str(subarea);
                end
                
                zmap = squeeze(PV(ian).zmap(nti,:,:,iv));
                zmapthresh = squeeze(PV(ian).clusterZmapThresh(nti, :,:,iv));
                contourf(PV(ian).time,PV(ian).frequency,zmap,40,'linecolor','none')
                hold on
                 [~,h] = contour(PV(ian).time,PV(ian).frequency,logical(zmapthresh),1);          
                  h.LineColor = 'black';
                set(gca,'ydir','normal','yscale','log');
                caxis(sf, 'auto')
                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',Pp.FontS,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end            
                            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('powerVarCorr %s %s', animal, PV(ian).expvars{iv});
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', Pp.FontS);
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                pconf = animaldef('Demetris');
                save_figure(sprintf('%s/powerVarCorr/',pconf{4}), 'powerVarCorr', sprtit)
                close all
            end
            close all;
        end
    end   
end