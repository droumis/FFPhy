
% need to get a value for the post .5:1 and pre -1:-.5 in theta band?
% whisker plot the post vs pre mean 
% then p

% plot mean per day
% plot mean per speed


pconf = paramconfig;
me = animaldef('Demetris');

run_ff = 0;
savedata = run_ff;
load_lfp = 0;
stack_lfp = 1;

load_stack = 0;
calc_AS = 1;
get_ripstate = 1;
load_ripstate = 0;
calc_mean_PWR = 1;
calc_ITPC = 1;

% loading and plotting only use 1 animal at a time
load_raw_pwr = 0; % unnormalized power (ntrode x sample x ripple x frequency)
plot_raw_pwr = 0;

load_mean_pwr = 0; % ripstate dB power means, baseline normed
plot_mean_pwr = 0;

load_raw_phase = 0;
plot_raw_phase = 0;

load_mean_itpc = 0;
plot_mean_itpc = 0;

pausefigs = 0;
savefigs = 1;

Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp.add_params = {'wtrackdays', 'excludeNoise','excludePriorFirstWell', '<4cm/s', ...
    'wavelets4-300Hz'};

Fp = load_filter_params(Fp, 'add_params', Fp.add_params);

Fp.uselfptype = 'eeg';
Fp.useripstates = {'onlywdays','rewarded', 'unrewarded', 'inbound' , 'outbound'};%, ...
% 'rewarded_inbound', 'unrewarded_inbound', 'rewarded_outbound', 'unrewarded_outbound'};
plot_frex = [8 25 100 200];
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
        sprintf('_%s%s', Fp.epochEnvironment))
end
%% load lfp
if load_lfp
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s%s', Fp.epochEnvironment));
end
%% stack LFP (ntrode x time x ripple)
if stack_lfp
    lfpstack = stack_riptriglfp(F);
    clear F
    save_data(lfpstack, Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack')
end
%% load lfpstack
if load_stack
    lfpstack = load_data(Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack', Fp.animals);
end

%% Analytic Signal
if calc_AS
    computeAnalyticSignal(lfpstack, 'waveSet', Fp.waveSet, 'saveOutput', 1, 'uselfptype', Fp.uselfptype);
end
%% RIPSTATE
if get_ripstate
    statefilters =  {'onlywdays', 'rewarded', 'unrewarded', 'inbound' , 'outbound'};
    ripstate = getStateFilters(lfpstack, 'statefilters', statefilters);
    save_data(ripstate, [pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack');
end
if load_ripstate
    ripstate = load_data([pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack', Fp.animals);
end
%% Compute and Save MEAN Power, ITPC
if calc_mean_PWR
    getPower(ripstate, Fp, 'uselfptype', Fp.uselfptype, 'ripstatetypes', ...
        Fp.useripstates, 'savepower', 1);
end
if calc_ITPC
    getITPC(ripstate, Fp, 'uselfptype', Fp.uselfptype, 'ripstatetypes', ...
        Fp.useripstates, 'saveresult', 1);
end
%% load power
if load_mean_pwr
    savedir = sprintf('%s/power/', pconf.andef{2});
    savestr=sprintf('/power_waveSet-%s_%s_%s',Fp.waveSet,Fp.uselfptype,Fp.epochEnvironment);
    pwr = load_data(savedir, savestr, Fp.animals);
end
%% load itpc
if load_mean_itpc
    savedir = sprintf('%s/itpc/', pconf.andef{2});
    savestr=sprintf('/itpc_waveSet-%s_%s_%s',Fp.waveSet,Fp.uselfptype,Fp.epochEnvironment);
    itpc = load_data(savedir, savestr, Fp.animals);
end

%% load raw power
if load_raw_pwr
    rawpwr = struct;
    for ian = 1:length(Fp.animals)
        animal = Fp.animals{ian};
        wp = getWaveParams(Fp.waveSet);
        rawpwr = load_data(sprintf('%s/analyticSignal/', me{2}), ...
            sprintf('AS_waveSet-%s_%s_power', wp.waveSet, Fp.uselfptype), animal);
    end
end
%% load raw itpc
if load_raw_phase
    rawphase = struct;
    for ian = 1:length(Fp.animals)
        animal = Fp.animals{ian};
        wp = getWaveParams(Fp.waveSet);
        rawphase = load_data(sprintf('%s/analyticSignal/', me{2}), ...
            sprintf('AS_waveSet-%s_%s_phase', wp.waveSet, Fp.uselfptype), animal);
    end
end
%%
if plot_raw_pwr
    Pp = load_plotting_params({'defaults', 'riptriglfp_perstatefrex_allntrodes'});
    for ian = 1:numel(Fp.animals)
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        anidx = find(strcmp({rawpwr.animal}, animal));
        ntrodes = lfpstack(anidx).ntrodes;
        dayep = [lfpstack(anidx).day lfpstack(anidx).epoch];
        
        % exclude invalid tets
        invalidtets = evaluatefilter(ntinfo, 'isequal($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        
        for fx = 1:length(plot_frex)
            frexidx = knnsearch(rawpwr.wp.frex', plot_frex(fx));
            knnfrex = round(rawpwr.wp.frex(frexidx));
            % plot per ripstate
            for rs = 1:length(Fp.useripstates)
                %         use_filts = find(any(cell2mat(cellfun(@(x) strcmp(x, lfpstack(anidx).filterfields), ...
                %             use_filters, 'un', 0)), 2));
                istate = Fp.useripstates{rs};
                istateidx = find(strcmp(istate, ripstate.statesetsfields));
                include_rips = ripstate.statesets(:,istateidx);
                %% ---- init fig----
                if savefigs && ~pausefigs
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                
                for nti = 1:length(ntrodes)
                    ntrode = ntrodes(nti);
                    if ismember(ntrode, invalidtets)
                        continue
                    end
                    
                    sf = subaxis(2,ceil(max(ntrodes)/2), nti, 'SpacingVert', Pp.SpVt, ...
                        'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                        Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    
                    exdayep = dayep(find(include_rips),:);
                    de = unique(exdayep, 'rows');
                    daybounds = find(diff(exdayep(:,1)));
                    epbounds = find(abs(diff(exdayep(:,2))));
%                     [d,i] = max(squeeze(double(pwr(anidx).pwr(11, :, find(include_rips),frexidx))));
%                     [j,k] = max(d);
%                     include_rips(k) = 0;
                    excld_stack = squeeze(double(rawpwr(anidx).pwr(nti, :, find(include_rips),frexidx)))';
                    idata2plot = trim2win(excld_stack, Fp.srate, Pp.pwin, ...
                        'dsamp', rawpwr(anidx).wp.dsamp);
                    %mididx = ceil(size(excld_stack,2)/2); % right now assumes center is rip start
                    ptime = linspace(-Pp.pwin(1),Pp.pwin(2),size(idata2plot,2));
                    z = idata2plot;
                    %                     m = nanmean(idata2plot,2);
                    %                     s = nanstd(idata2plot, [], 2);
                    %                     z = (idata2plot-m)./s;
                    imagesc(ptime, 1:size(z,1), z)
                    colormap(parula)
                    caxis(sf, 'auto')
                    
                    line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color',[.9 .9 .9])
                    line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
                    line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
                    title(sprintf('%d', ntrode), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                        'FontName', Pp.FontNm)
                    %                         caxis([-1 1])
                    xlabel('time s', 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', ...
                        Pp.FontNm)
                    ylabel('ripnum (day-b epoch-w)','FontSize',Pp.FontS, ...
                        'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                    if nti ~= 1
                        xlabel('')
                        ylabel('')
                        set(gca, 'ytick', []);
                        set(gca, 'xtick', []);
                    end
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('%s %s %s %dHz %s power', animal, istate, Fp.uselfptype, ...
                    knnfrex, Fp.epochEnvironment);
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                spylabel = text(.02, .5, sprintf('ripnum'), 'Parent', sprtitleax, 'Units', ...
                    'normalized', 'Rotation', 90);
                set(spylabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center', 'FontSize', 12);
                
                spxlabel = text(.5, .02, sprintf('time'), 'Parent', sprtitleax, 'Units', ...
                    'normalized');
                set(spxlabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center', 'FontSize', 12);
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    pconf = animaldef('Demetris');
                    save_figure(sprintf('%s/wavepower/',pconf{4}), 'wavepower', sprtit)
                    close all
                end
            end
        end
    end
end

%% plot power
if plot_mean_pwr
    Pp = load_plotting_params({'defaults', 'power'});
    wp = getWaveParams('4-300Hz');
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo, 'area');
        matidx = unique(den.index(:,3));
        anidx = find(strcmp({pwr.animal}, animal));
        for co = 1:length(Fp.useripstates)
            if savefigs && ~pausefigs
                close all
                ifig =figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            numcols = 5;
            numrows = ceil(length(ntrodes) / numcols);
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
                ntidx = find(matidx == nt);
                idata2plot = squeeze(pwr(anidx).meandbpower{co}.pwr_mean_db(ntidx,:,:))';
                idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                    'dsamp', pwr(anidx).wp.dsamp);
                time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                contourf(sf, time, wp.frex, idata2plot, Pp.contourRes, ...
                    'linecolor','none');
                set(gca,'ydir','normal','yscale','log');
                
                colormap(Pp.usecolormap)
                caxis(sf, 'auto')
                
                hold on
                % thresholded single pix zmask
                if ~isempty(fieldnames(pwr(anidx).meandbpower{co}.permt))
                    zmask2plot = squeeze(pwr(anidx).dbpower{co}.permt.threshmean(ntidx,:,:))'; % ntrode is in 3rd dim here
                    zmask2plot = trim2win(zmask2plot, Fp.srate, Pp.pwin, 'dsamp', dsamp);
                    [~,h] = contour(sf, time, wp.frex, logical(zmask2plot), 1);
                    h.LineColor = 'black';
                end
                hold on;
                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s %s %s %s norm mean dbPower', animal, 'allnts', ...
                Fp.uselfptype, Fp.useripstates{co});
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 16);
            
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

if plot_raw_phase
    Pp = load_plotting_params({'defaults', 'riptriglfp_perstatefrex_allntrodes'});
    for ian = 1:numel(Fp.animals)
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        anidx = find(strcmp({rawphase.animal}, animal));
        ntrodes = lfpstack(anidx).ntrodes;
        dayep = [lfpstack(anidx).day lfpstack(anidx).epoch];
        
        % exclude invalid tets
        invalidtets = evaluatefilter(ntinfo, 'isequal($valid, ''no'')');
        invalidtets = unique(invalidtets(:,3));
        
        for fx = 1:length(plot_frex)
            frexidx = knnsearch(rawphase.wp.frex', plot_frex(fx));
            knnfrex = round(rawphase.wp.frex(frexidx));
            % plot per ripstate
            for rs = 1:length(Fp.useripstates)
                %         use_filts = find(any(cell2mat(cellfun(@(x) strcmp(x, lfpstack(anidx).filterfields), ...
                %             use_filters, 'un', 0)), 2));
                istate = Fp.useripstates{rs};
                istateidx = find(strcmp(istate, ripstate.statesetsfields));
                include_rips = ripstate.statesets(:,istateidx);
                %% ---- init fig----
                if savefigs && ~pausefigs
                    close all
                    ifig = figure('Visible','off','units','normalized','position', ...
                        Pp.position);
                else
                    ifig = figure('units','normalized','position',Pp.position);
                end
                set(gcf,'color','white')
                
                for nti = 1:length(ntrodes)
                    ntrode = ntrodes(nti);
                    if ismember(ntrode, invalidtets)
                        continue
                    end
                    
                    sf = subaxis(2,ceil(max(ntrodes)/2), nti, 'SpacingVert', Pp.SpVt, ...
                        'SpacingHoriz', Pp.SpHz, 'MarginLeft', Pp.MgLt, 'MarginRight', ...
                        Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
                    
                    exdayep = dayep(find(include_rips),:);
                    de = unique(exdayep, 'rows');
                    daybounds = find(diff(exdayep(:,1)));
                    epbounds = find(abs(diff(exdayep(:,2))));
%                     [d,i] = max(squeeze(double(pwr(anidx).pwr(11, :, find(include_rips),frexidx))));
%                     [j,k] = max(d);
%                     include_rips(k) = 0;
                    excld_stack = squeeze(double(rawphase(anidx).ph(nti, :, find(include_rips),frexidx)))';
                    idata2plot = trim2win(excld_stack, Fp.srate, Pp.pwin, ...
                        'dsamp', rawphase(anidx).wp.dsamp);
                    %mididx = ceil(size(excld_stack,2)/2); % right now assumes center is rip start
                    ptime = linspace(-Pp.pwin(1),Pp.pwin(2),size(idata2plot,2));
                    z = idata2plot;
                    %                     m = nanmean(idata2plot,2);
                    %                     s = nanstd(idata2plot, [], 2);
                    %                     z = (idata2plot-m)./s;
                    imagesc(ptime, 1:size(z,1), z)
                    colormap(parula)
                    caxis(sf, 'auto')
                    
                    line([-Pp.pwin(1) Pp.pwin(2)], [epbounds'; epbounds'], 'color',[.9 .9 .9])
                    line([-Pp.pwin(1) Pp.pwin(2)], [daybounds'; daybounds'], 'color', [0 0 0])
                    line([0 0], [1 size(z,1)], 'color', [0 0 0], 'linestyle', '--')
                    title(sprintf('%d', ntrode), 'FontSize',Pp.FontS, 'FontWeight',Pp.FontW, ...
                        'FontName', Pp.FontNm)
                    %                         caxis([-1 1])
                    xlabel('time s', 'FontSize',Pp.FontS,'FontWeight',Pp.FontW,'FontName', ...
                        Pp.FontNm)
                    ylabel('ripnum (day-b epoch-w)','FontSize',Pp.FontS, ...
                        'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
                    if nti ~= 1
                        xlabel('')
                        ylabel('')
                        set(gca, 'ytick', []);
                        set(gca, 'xtick', []);
                    end
                end
                %% super
                sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
                sprtit = sprintf('%s %s %s %dHz %s phase', animal, istate, Fp.uselfptype, ...
                    knnfrex, Fp.epochEnvironment);
                iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
                set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center','FontSize', 12);
                spylabel = text(.02, .5, sprintf('ripnum'), 'Parent', sprtitleax, 'Units', ...
                    'normalized', 'Rotation', 90);
                set(spylabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center', 'FontSize', 12);
                
                spxlabel = text(.5, .02, sprintf('time'), 'Parent', sprtitleax, 'Units', ...
                    'normalized');
                set(spxlabel,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                    'horizontalAlignment', 'center', 'FontSize', 12);
                %% ---- pause, save figs ----
                if pausefigs
                    pause
                end
                if savefigs
                    pconf = animaldef('Demetris');
                    save_figure(sprintf('%s/itpc/',pconf{4}), 'itpc', sprtit)
                    close all
                end
            end
        end
    end
end
%% plot ITPC
if plot_mean_itpc
    Pp = load_plotting_params({'defaults', 'power'});
    wp = getWaveParams('4-300Hz');
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        aninfo = animaldef(animal);
        ntinfo = loaddatastruct(aninfo{2}, animal, 'tetinfo');
        ntrodes = evaluatefilter(ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo, 'area');
        matidx = unique(den.index(:,3));
        anidx = find(strcmp({pwr.animal}, animal));
        for co = 1:length(Fp.useripstates)
            if savefigs && ~pausefigs
                close all
                ifig =figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            numcols = 5;
            numrows = ceil(length(ntrodes) / numcols);
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
                ntidx = find(matidx == nt);
                idata2plot = squeeze(itpc(anidx).ITPC{co}.ITPC(ntidx,:,:))';
                idata2plot = trim2win(idata2plot, Fp.srate, Pp.pwin, ...
                    'dsamp', itpc(anidx).wp.dsamp);
                time = linspace(-Pp.pwin(1), Pp.pwin(2), length(idata2plot(1,:)));
                contourf(sf, time, wp.frex, idata2plot, Pp.contourRes, ...
                    'linecolor','none');
                set(gca,'ydir','normal','yscale','log');
                colormap(Pp.usecolormap)
                caxis(sf, 'auto')
                
                hold on
                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 8)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', Pp.FontNm)
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s %s %s %s itpc', animal, 'allnts', ...
                Fp.uselfptype, Fp.useripstates{co});
            iStitle = text(.5, .98, {sprtit}, 'Parent', sprtitleax, 'Units', 'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center','FontSize', 16);
            
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