

% combine the mec ntrodes into areas and re-plot the mean/itpc per animal / area
% - post vs pre
% - across days


load_pwr = 0;
plot_pwr = 1;

load_itpc = 0;
plot_itpc = 1;

plot_prepost = 0;

pausefigs = 0;
savefigs = 1;

Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'};
Fp.add_params = {'wtrack', 'wavelets4-300Hz', 'excludeNoise','excludePriorFirstWell'};
Fp.filtfunction = 'dfa_riptriglfp';
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
pconf = paramconfig;

Fp.uselfptype = 'eeggnd';
Fp.useripstates = {'all'};
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
