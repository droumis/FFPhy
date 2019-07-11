% from loading the lfp to printing a power figure is 20 min / animal
% 2.5 hours for all 7 animals if run on the same computer
% if i run 3 animals on virga01 and 4 on typhoon, only 1.5 hours total

load_stack = 1;
calc_AS = 1;

get_ripstate = 0;
load_ripstate = 1;
calc_PWR = 1;
calc_ITPC = 1;
run_permtest = 0;

load_pwr = 0;
plot_pwr = 1;

load_itpc = 0;
plot_itpc = 1;

plot_prepost = 0;

pausefigs = 0;
savefigs = 1;

Fp.animals = {'JZ2', 'JZ3', 'JZ4'};
Fp.add_params = {'wtrack', 'wavelets4-300Hz', 'excludeNoise','excludePriorFirstWell'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
pconf = paramconfig;

%% load lfpstack ~ 1 minute on
if load_stack
    lfpstack = load_data(Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack', Fp.animals);
end
%% Calculate and Save Analytic Signal ~ 15 minutes / animal on virga01
% virga01. D10 compute AS: 2 minutes, saving phase,power each: 6 minutes (12Gb)
% 437 seconds to make cmwfft on typhoon. 275 seconds to compute as, 1084 s
% to save each
if calc_AS
    computeAnalyticSignal(lfpstack, 'waveSet', Fp.waveSet, 'overwrite', 1, ...
        'saveAnalyticSignal', 1, 'uselfptype', Fp.uselfptype);
end
%% get behavioral state for each rip
if get_ripstate
    ripstate = getStateFilters(lfpstack);
    save_data(ripstate, [pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack');
end
if load_ripstate
    ripstate = load_data([pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack', Fp.animals);
end
%% Compute and Save Power ~ 2 minutes on virga01
if calc_PWR
    pwr = getPower(ripstate, Fp, 'uselfptype', Fp.uselfptype, 'ripstatetypes', ...
        Fp.useripstates, 'run_permutation_test', run_permtest, 'savepower', 1);
end
if calc_ITPC
    itpc = getITPC(ripstate, Fp, 'uselfptype', Fp.uselfptype, 'ripstatetypes', ...
        Fp.useripstates, 'run_permutation_test', run_permtest, 'saveresult', 1);
end
%% load power
if load_pwr
    savedir = sprintf('%s/power/', pconf.andef{2});
    savestr=sprintf('/power_waveSet-%s_%s_%s',Fp.waveSet,Fp.uselfptype,Fp.epochEnvironment);
    pwr = load_data(savedir, savestr, Fp.animals);
end
%% load itpc
if load_itpc
    savedir = sprintf('%s/itpc/', pconf.andef{2});
    savestr=sprintf('/itpc_waveSet-%s_%s_%s',Fp.waveSet,Fp.uselfptype,Fp.epochEnvironment);
    itpc = load_data(savedir, savestr, Fp.animals);
end

%% plot power
if plot_pwr
    ntinfo = struct;
    for ian = 1:length(Fp.animals)
        ntinfo(ian).animal = Fp.animals{ian};
        andef = animaldef(ntinfo(ian).animal);
        ntinfo(ian).ntinfo = loaddatastruct(andef{2}, ntinfo(ian).animal, 'tetinfo');
    end
    Pp = load_plotting_params({'defaults', 'power'});
    wp = getWaveParams('4-300Hz');
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        anidx = find(strcmp({pwr.animal}, animal));
        invalidtets = evaluatefilter(ntinfo(anidx).ntinfo, 'strcmp($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        areas = {'ca1', 'mec', 'ref'};
        den = cellfetch(ntinfo(anidx).ntinfo, '');
        matidx = unique(den.index(:,3));
        for ar = 1:length(areas) % for each area
            areantrodes = evaluatefilter(ntinfo(anidx).ntinfo, ...
                sprintf('isequal($area, ''%s'')', areas{ar}));
            areantrodes = unique(areantrodes(~ismember(areantrodes(:,3), invalidtets),3));
            % for each condition state
            for co = 1:length(Fp.useripstates)
                if savefigs && ~pausefigs
                    close all
                    ifig =figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                for nti = 1:length(areantrodes)
                    sf = subaxis(4,5,nti, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                        Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    nt = areantrodes(nti);
                    ntidx = find(matidx == nt);
                    idata2plot = squeeze(pwr(anidx).meandbpower{co}.pwr_mean_db(ntidx,:,:))';
                    idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                        'dsamp', pwr(anidx).wp.dsamp);
                    time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                    contourf(sf, time, wp.frex, idata2plot, Pp.contourRes, ...
                        'linecolor','none');
                    set(gca,'ydir','normal','yscale','log');
                    
                    if strcmp(areas{ar}, 'mec')
                        coloraxis = [-1 1]; %[-5 50];
                    elseif strcmp(areas{ar}, 'ca1')
                        coloraxis = [-2 5]; %[-5 50];
                    end
                    caxis(coloraxis)
                    colormap(Pp.usecolormap)
                    hold on
                    % thresholded single pix zmask
                    if ~isempty(fieldnames(pwr(anidx).meandbpower{co}.permt))
                        zmask2plot = squeeze(pwr(anidx).dbpower{co}.permt.threshmean(ntidx,:,:))'; % ntrode is in 3rd dim here
                        zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', dsamp);
                        [~,h] = contour(sf, time, wp.frex, logical(zmask2plot), 1);
                        h.LineColor = 'black';
                    end
                    hold on;
                    subarea = ntinfo(anidx).ntinfo{1}{1}{nt}.subarea;
                    if isnumeric(subarea)
                        subarea = num2str(subarea);
                    end
                    title(sprintf('%s %s, nt%d',areas{ar},subarea,nt), 'FontSize',Pp.FontS,...
                        'FontWeight',Pp.FontW, 'FontName', ...
                        Pp.FontNm)
                    if mod(nti+4,5) == 0
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
                sprtit = sprintf('%s %s mean dbPower p01zmap %s %s %s %s', animal, areas{ar}, ...
                    Fp.add_params{1}, Fp.add_params{2});
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    pconf = animaldef('Demetris');
                    save_figure(sprintf('%s/wavepower/',pconf{4}), 'wavepower', sprtit)
                    close all
                end
                close all;
                %                 end
            end
        end
    end
end

%% Plot ITPC
if plot_itpc
    ntinfo = struct;
    for ian = 1:length(Fp.animals)
        ntinfo(ian).animal = Fp.animals{ian};
        andef = animaldef(ntinfo(ian).animal);
        ntinfo(ian).ntinfo = loaddatastruct(andef{2}, ntinfo(ian).animal, 'tetinfo');
    end
    Pp = load_plotting_params({'defaults', 'power'});
    wp = getWaveParams('4-300Hz');
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        anidx = find(strcmp({itpc.animal}, animal));
        invalidtets = evaluatefilter(ntinfo(anidx).ntinfo, 'strcmp($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        areas = {'ca1', 'mec', 'ref'};
        den = cellfetch(ntinfo(anidx).ntinfo, '');
        matidx = unique(den.index(:,3));
        for ar = 1:length(areas) % for each area
            areantrodes = evaluatefilter(ntinfo(anidx).ntinfo, ...
                sprintf('isequal($area, ''%s'')', areas{ar}));
            areantrodes = unique(areantrodes(~ismember(areantrodes(:,3), invalidtets),3));
            for co = 1:length(Fp.useripstates)
                if savefigs && ~pausefigs
                    close all
                    ifig =figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                for nti = 1:length(areantrodes)
                    sf = subaxis(4,5,nti, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                        Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    nt = areantrodes(nti);
                    ntidx = find(matidx == nt);
                    idata2plot = squeeze(itpc(anidx).ITPC{co}.ITPC(ntidx,:,:))';
                    idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                        'dsamp', itpc(anidx).wp.dsamp);
                    time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                    contourf(sf, time, wp.frex, idata2plot, Pp.contourRes, ...
                        'linecolor','none');
                    set(gca,'ydir','normal','yscale','log');
                    caxis('auto')
                    colormap(Pp.usecolormap)
                    hold on
                    % thresholded single pix zmask
                    if ~isempty(fieldnames(itpc(anidx).ITPC{co}.permt))
                        zmask2plot = squeeze(itpc(anidx).ITPC{co}.permt.threshmean(ntidx,:,:))'; % ntrode is in 3rd dim here
                        zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', dsamp);
                        [~,h] = contour(sf, time, wp.frex, logical(zmask2plot), 1);
                        h.LineColor = 'black';
                    end
                    hold on;
                    subarea = ntinfo(anidx).ntinfo{1}{1}{nt}.subarea;
                    if isnumeric(subarea)
                        subarea = num2str(subarea);
                    end
                    title(sprintf('%s %s, nt%d',areas{ar},subarea,nt), 'FontSize',Pp.FontS,...
                        'FontWeight',Pp.FontW, 'FontName', ...
                        Pp.FontNm)
                    if mod(nti+4,5) == 0
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
                sprtit = sprintf('%s %s ITPC %s %s %s %s', animal, areas{ar}, ...
                    Fp.add_params{1}, Fp.add_params{2});
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    pconf = animaldef('Demetris');
                    save_figure(sprintf('%s/itpc/',pconf{4}), 'itpc', sprtit)
                    close all
                end
                close all;
            end
        end
    end
end

%% plot_prepost
if plot_prepost
    ntinfo = struct;
    for ian = 1:length(Fp.animals)
        ntinfo(ian).animal = Fp.animals{ian};
        andef = animaldef(ntinfo(ian).animal);
        ntinfo(ian).ntinfo = loaddatastruct(andef{2}, ntinfo(ian).animal, 'tetinfo');
    end
    Pp = load_plotting_params({'defaults', 'power'});
    wp = getWaveParams('4-300Hz');
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        anidx = find(strcmp({pwr.animal}, animal));
        invalidtets = evaluatefilter(ntinfo(anidx).ntinfo, 'strcmp($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        areas = {'ca1', 'mec', 'ref'};
        den = cellfetch(ntinfo(anidx).ntinfo, '');
        matidx = unique(den.index(:,3));
        for ar = 1:length(areas) % for each area
            areantrodes = evaluatefilter(ntinfo(anidx).ntinfo, ...
                sprintf('isequal($area, ''%s'')', areas{ar}));
            areantrodes = unique(areantrodes(~ismember(areantrodes(:,3), invalidtets),3));
            for co = 1:length(Fp.useripstates)
                if savefigs && ~pausefigs
                    close all
                    ifig =figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                for nti = 1:length(areantrodes)
                    sf = subaxis(4,5,nti, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                        'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', ...
                        Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    nt = areantrodes(nti);
                    ntidx = find(matidx == nt);
                    predata = squeeze(pwr(anidx).meandbpower{co}.premean(ntidx,:,:));
                    postdata = squeeze(pwr(anidx).meandbpower{co}.postmean(ntidx,:,:));
                    plot(mean(postdata-predata))
                    set(gca,'xtick', 1:3:size(postdata,2))
                    set(gca,'XTickLabel',round(wp.frex(1:3:end)))
                    xlabel('frequency')
                    ylabel('mean post-pre win')
                    subplot(2,2,1)
                    [a,i]=max(postdata-predata);
                    histogram(); %, [], 1:length(postdata))
                    axis tight
                    subplot(2,2,2)
                    %                     boxplot(
                    %                     allData = {Data1; Data2; Data3};
                    %                     h = boxplot([allData{:}],group);
                    %                     set(h, 'linewidth' ,2)
                    %                     set(gca,'XTickLabel', {'Data1'; 'Data2'; 'Data3'})
                    
                    title(sprintf('%s %s, nt%d',areas{ar},subarea,nt), 'FontSize',Pp.FontS,...
                        'FontWeight',Pp.FontW, 'FontName', ...
                        Pp.FontNm)
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('%s %s mean dbPower p01zmap %s %s %s %s', animal, areas{ar}, ...
                    Fp.add_params{1}, Fp.add_params{2});
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    pconf = animaldef('Demetris');
                    save_figure(sprintf('%s/wavepower/',pconf{4}), 'wavepower', sprtit)
                    close all
                end
                close all;
                %                 end
            end
        end
    end
end

