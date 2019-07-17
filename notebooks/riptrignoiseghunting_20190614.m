

% the goal of all this is to get a per area (mec sup, mec deep, ca1) riptrig power and itpc plot
% this means combined across animals and days and ntrodes

% need to detect and remove the noisy ripples from all the animals
% the noise is apparent in the wavelet power and itpc plots
% do this with unreferenced lfp since, the referencing is problematic for a
% couple of the animals still
% the referencing is probably also something i need to keep debugging.. 

% load the power stack. plot time vs rip # vs power for the nearest
% frequency to the band where the noise is most apparent in the power mean
% plots

% then observe which are the problematic noisy ripples.. maybe by sorting
% by the max in that specific time-frequency bin that seems most indicative
% of the issue
pconf = paramconfig;

load_stack = 0;
calc_AS = 0;

get_ripstate = 1;
load_ripstate = 0;
calc_PWR = 1;
calc_ITPC = 1;
run_permtest = 0;

load_pwr = 1;
plot_pwr = 1;

load_itpc = 1;
plot_itpc = 1;

plot_prepost = 0;

pausefigs = 0;
savefigs = 1;

Fp.animals = {'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp.add_params = {'wtrackdays', 'excludeNoise','excludePriorFirstWell', '<4cm/s', ...
    'wavelets4-300Hz'}; %'correcttrials', 

Fp = load_filter_params(Fp, 'add_params', Fp.add_params);

Fp.uselfptype = 'eeggnd';
Fp.useripstates = {'all', 'onlywdays', 'rewarded', 'unrewarded', 'inbound' , 'outbound', ...
    'rewarded_inbound', 'unrewarded_inbound', 'rewarded_outbound', 'unrewarded_outbound'};
%% extract from rec
%% adjust timestamps
%% reference lfp
%% run ff get riptrig lfp 
%% stack ff result per animal

%% load lfpstack
if load_stack
    lfpstack = load_data(Fp.paths.resultsDirectory, 'riptriglfpstack_wtrack', Fp.animals);
end
%% get behavioral state for each rip
if get_ripstate
    ripstate = getStateFilters(lfpstack);
    save_data(ripstate, [pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack');
end
if load_ripstate
    ripstate = load_data([pconf.andef{2}, 'ripstate/'], 'ripstate_wtrack', Fp.animals);
end

%% Compute and Save Power
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

for ian = 1:length(Fp.animals)
    ntinfo(ian).animal = Fp.animals{ian};
    andef = animaldef(ntinfo(ian).animal);
    ntinfo(ian).ntinfo = loaddatastruct(andef{2}, ntinfo(ian).animal, 'tetinfo');
end
% RUN again because i doubled the baseline period in hopes of fixing the weirdness with JZ2

% so i really want to be able to identify which ripples are noisy.. and
% exclude them in a principled way.. so i need to be able to see across
% ripples in the freqeucny band where the noise is very apparent.. 
% for D12 that's all the high frequencies.. so just use the highest one.. i
% think 300. 
% once i plot the noise freq specific heatrast for all the ntrodes and get a sense of which
% ripples are noisy.. i want to identify which ones, day, ep, timestamp.. 
% i could maybe sort by the highest val within that noise specific time
% window and then look at the individual ripples at the top in the per rip
% plots to see if there's zsomething weird.. 

% i want to plot a sort of ripple condition/state key that shows which
% ripples were included.. this could use a rip heatrast stack from one
% ntrode.. and/or plotted along with the posperform plots.. although i'd
% need to condense that a lot somehow.. i think i should plot these
% seperately and just make sure i label them with what condition they are.

% having this key alongside the power/itpc mean results will help me build
% confidence about the results i'm seeing. 
% other important information would be the total number of ripples.. since
% we probably want this to be close for any two conditions we want to
% compare head to head.. 

% i also think that i want to plot the actual plots themselves in a layout
% that reflects their position in the cannula. but that might be for
% another time?

% in terms of 
%% plot power
if plot_pwr
    Pp = load_plotting_params({'defaults', 'power'});
    wp = getWaveParams('4-300Hz');
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        
        ntinfanidx = find(strcmp({ntinfo.animal}, animal));
        ntrodes = evaluatefilter(ntinfo(ntinfanidx).ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo(ntinfanidx).ntinfo, 'area');
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
                area = ntinfo(ntinfanidx).ntinfo{1}{1}{nt}.area;
                subarea = ntinfo(ntinfanidx).ntinfo{1}{1}{nt}.subarea;
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
%                 colorbar
                
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
%                 xlabel('time s', 'FontSize',8,'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
%                 ylabel('freq Hz','FontSize',8, 'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
%                 else
%                     set(gca, 'xlabel', [])
%                     set(gca, 'ylabel', [])
%                 end
%                                 
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s %s norm mean dbPower %s %s %s %s', animal, 'allnts', ...
                Fp.uselfptype, Fp.useripstates{co}, Fp.add_params{1}, Fp.add_params{2});
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
 
%% plot ITPC
if plot_itpc
    Pp = load_plotting_params({'defaults', 'power'});
    wp = getWaveParams('4-300Hz');
    for ian = 1:numel(Fp.animals) % for each animal
        animal = Fp.animals{ian};
        
        ntinfanidx = find(strcmp({ntinfo.animal}, animal));
        ntrodes = evaluatefilter(ntinfo(ntinfanidx).ntinfo, 'strcmp($valid, ''yes'')');
        ntrodes = unique(ntrodes(:,3));
        den = cellfetch(ntinfo(ntinfanidx).ntinfo, 'area');
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
                area = ntinfo(ntinfanidx).ntinfo{1}{1}{nt}.area;
                subarea = ntinfo(ntinfanidx).ntinfo{1}{1}{nt}.subarea;
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
%                 colorbar
                hold on                
                ytickskip = 2:4:wp.numfrex;
                set(gca,'ytick', round(wp.frex(ytickskip)), 'FontSize', 7)
                title(sprintf('%s%s nt%d',area,subarea,nt), 'FontSize',14,...
                    'FontWeight',Pp.FontW, 'FontName', ...
                    Pp.FontNm)
%                 xlabel('time s', 'FontSize',8,'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
%                 ylabel('freq Hz','FontSize',8, 'FontWeight',Pp.FontW,'FontName', Pp.FontNm)
%                 else 
%                     set(gca, 'xlabel', [])
%                     set(gca, 'ylabel', [])
%                 end
%                                 
                yl = ylim;
                line([0 0], yl, 'Color', [0.8 0.8 0.8],'LineStyle','--', 'LineWidth', 1);
            end
            %% super
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('%s %s itpc %s %s %s %s', animal, 'allnts', ...
                Fp.uselfptype, Fp.useripstates{co}, Fp.add_params{1}, Fp.add_params{2});
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
            %                 end
        end
    end
end
