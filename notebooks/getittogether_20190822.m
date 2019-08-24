
% get speed around ripple times

animal = 'D10';
andef = animaldef(animal);
swr = loaddatastruct(andef{2}, animal, 'ca1rippleskons');
pos = loaddatastruct(andef{2}, animal, 'pos');
%%
day = 6;
epoch = 2;
%%
velcol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'vel-loess'));
timecol = find(strcmp(strsplit(pos{day}{epoch}.fields,' '), 'time'));
timevel = pos{day}{epoch}.data(:, [timecol velcol]);
%% Velocity CDF
subplot(2,2,1)
histogram(timevel(:,2), 'Normalization', 'cdf', ...
    'FaceColor', 'black', 'LineStyle', 'none');
title('Velocity CDF')
ylabel('Probability')
xlabel('Log10 Velocity cm/s')
set(gca, 'xscale', 'log')
set(gca, 'xtick', [0 2 4 6 8 10 20 40])
set(gca, 'ytick', 0:.2:1)
% axis tight
xlim([0 60])
%% Inter-SWR starttime Interval CDF..
subplot(2,2,2)
histogram(diff(swr{day}{epoch}{1}.starttime), 100, 'Normalization', 'cdf',...
    'FaceColor', [0 0 .3], 'LineStyle', 'none');
title('Inter-SWR Interval CDF')
xlabel('Log10 time s')
ylabel('Probability')
set(gca, 'xscale', 'log')
set(gca, 'xtick', [0 2 4 6 8 10 20 40])
set(gca, 'ytick', 0:.2:1)
% axis tight
xlim([0 60])
ylim([0 1])
%% SWR Duration CDF
subplot(2,2,3)
% swrstr.data = diff(swr{day}{epoch}{1}.starttime);
% swrstr.descript = 'swrstarttime';
% ac = acorr(swrstr, 10, 200);
% ac.data(ac.data(:,1)==0,:) = []; % drop 0-lag
% b = bar(ac.data(:,1)/1000, ac.data(:,2));
% ac = acf(round(diff(swrstr.data)*1000), 100);
% xlabel()
swrdur = swr{day}{epoch}{1}.endtime - swr{day}{epoch}{1}.starttime;
histogram(swrdur*1000, 100, 'Normalization', 'cdf',...
    'FaceColor', [0 .3 .3], 'LineStyle', 'none');
title('SWR Duration CDF')
xlabel('Log10 time ms')
ylabel('Probability')
set(gca, 'xscale', 'log')
set(gca, 'yscale', 'linear')
set(gca, 'xtick', [0 50 100 200])
axis tight
% set(gca, 'ytick', 0:.2:1)

%% Peri swr speed histogram
subplot(2,2,4)
%% ripLFP [ntrode time rip]
Fp.animals = {'D10', 'D12', 'D13', 'JZ1', 'JZ2', 'JZ3', 'JZ4'}; %, 'JZ4'}; %;{'D10'}; %
% Fp.animals = {'JZ3'};
Fp.filtfunction = 'dfa_getPeriEventVelocity';
Fp.add_params = {'wtrackdays','excludeNoise','excludePriorFirstWell','<4cm/s', ...
    'excludeAfterLastWell'};
Fp = load_filter_params(Fp, 'add_params', Fp.add_params);
pconf = paramconfig;
%%
make_swrLFP = 0;
save_swrLFP = 0;
load_swrLFP = 1;
if make_swrLFP
    F = createfilter('animal', Fp.animals, 'epochs', Fp.epochfilter, ...
        'excludetime', Fp.timefilter, 'iterator', Fp.iterator);
    F = setfilterfunction(F, Fp.filtfunction, Fp.datatypes, Fp.options{:});
    F = runfilter(F); for d = 1:length(F); F(d).datafilter_params = Fp; end; end
if save_swrLFP
    save_data(F, Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, 'filetail',...
        sprintf('_%s', Fp.epochEnvironment)); end
if load_swrLFP
    F = load_data(Fp.paths.filtOutputDirectory, Fp.paths.filenamesave, ...
        Fp.animals, 'filetail', sprintf('_%s', Fp.epochEnvironment)); end
% %%
% allrips = cell2mat(permute({F.output{1}.data},[3 1 2]));
% velmean = squeeze(median(allrips(9,:,:),3));
% velsem = squeeze(std(allrips(9,:,:),[],3))/sqrt(size(allrips,3));
%%
% subplot(2,2,4)
% shadedErrorBar(F.output{1}(1).time, velmean, velsem)
% title('per-SWR velocity')
% xlabel('time s')
% ylabel('velocity')
% set(gca, 'xscale', 'linear')
% set(gca, 'yscale', 'linear')
% % set(gca, 'xtick', [0 50 100 200])
% yl = ylim;
% % ylim([0,yl(2)])
% axis tight
%%
% - GOOD.. now do all the animals and make a per-SWR vel plot per animal
% - then make a few plots of immobility time (with speed, pos) with MEC LFP
% (sup, deep labled).. 
plotout = 1;
pausefigs = 0;
savefigs = 1;
if plotout; Pp=load_plotting_params({'defaults','riptrigvel'});
    if savefigs && ~pausefigs; close all
        ifig = figure('Visible','off','units','normalized','position', Pp.position);
    else; ifig = figure('units','normalized','position',Pp.position); end
    set(gcf,'color','white')
        for ani = 1:length(F)
            sf = subaxis(2,7,ani, 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
                    'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
            allrips = cell2mat(permute({F(ani).output{1}.data},[3 1 2]));
            velmean = squeeze(mean(allrips(9,:,:),3));
            velstd = squeeze(std(allrips(9,:,:),[],3)); %
            velsem = velstd/sqrt(size(allrips,3));
%             shadedErrorBar(F(ani).output{1}(1).time, velmean, velstd)
            hold on;
            shadedErrorBar(F(ani).output{1}(1).time, velmean, velsem);%, {'color', [0 0 .2 .2]})
            title(F(ani).animal{3});
            if ani == 1
                xlabel('time s')
                ylabel('velocity mean, sem')
            end
            axis tight
            yl = ylim;
            ylim([0 yl(2)]);
            xl = xlim;
            line(xl, [4 4], 'color', [0 0 .5], 'linewidth', 2)
            x = [-1 -.5 -.5 -1];
            y = [yl(1) yl(1) yl(2) yl(2)];
            patch(x,y,'black', 'FaceAlpha', .1, 'LineStyle','none')
            x = [.5 1 1 .5];
            y = [yl(1) yl(1) yl(2) yl(2)];
            patch(x,y,'black', 'FaceAlpha', .1,'LineStyle','none')
            axis tight
            
            % plot mean/std bars of pre (-1:-.5) and post (.5:1) and do ttest
            sf = subaxis(2,7,ani+length(F), 'SpacingVert', Pp.SpVt, 'SpacingHoriz', Pp.SpHz, ...
            'MarginLeft', Pp.MgLt, 'MarginRight', Pp.MgRt, 'MarginTop', Pp.MgTp, 'MarginBottom', Pp.MgBm);
            t = knnsearch(F(ani).output{1}(1).time', [-1 -.5 .5 1]');    
            preM = squeeze(mean(mean(allrips(9,t(1):t(2),:),2),3));
            postM = squeeze(mean(mean(allrips(9,t(3):t(4),:),2),3));
            preS = squeeze(std(mean(allrips(9,t(1):t(2),:),2),[],3));
            postS = squeeze(std(mean(allrips(9,t(3):t(4),:),2),[],3));
%             a = reshape(allrips(9,t(1):t(2),:),[],1);
%             b = reshape(allrips(9,t(3):t(4),:),[],1);
            
            bar([1 2], [preM postM], 'facecolor', [.2 .2 .2])
            hold on
            snames={'pre'; 'post'};
            set(gca,'xticklabel',snames)
            er = errorbar([preM postM], [preS postS]);
            er.Color = [0 0 0];                         
            er.LineStyle = 'none';
            hold off
            if ani == 1
                ylabel('velocity mean, std')
            end
        end
        % ---- super axis -----
        sprax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
        sprtit = sprintf('Peri-SWR Velocity %s', Fp.epochEnvironment);
        iStitle = text(.5, .98, {sprtit}, 'Parent', sprax, 'Units', 'normalized');
        set(iStitle,'FontWeight','bold','FontName','Arial','horizontalAlignment','center','FontSize',12);
        % ---- pause, save figs ----
        if pausefigs; pause; end
        if savefigs; save_figure([pconf.andef{4},'/' Fp.filtfunction '/'],sprtit); end
end
%%

% swrstr.data = diff(swr{day}{epoch}{1}.starttime);
% swrstr.descript = 'swrstarttime';
% ac = acorr(swrstr, 10, 200);
% ac.data(ac.data(:,1)==0,:) = []; % drop 0-lag
% b = histogram((:,1)/1000, ac.data(:,2));
% axis tight




%%
      x = randn(1000,1);         % 1000 Gaussian deviates ~ N(0,1)
      y = filter([1 -1 1],1,x);  % Create an MA(2) process
      [acf,lags,bounds] = autocorr(y);
      stem(lags,acf); xlabel('Lag'); ylabel('\rho(k)');
      hold on;
      h = line(lags,bounds(1)*ones(length(acf),1));
      h1 = line(lags,bounds(2)*ones(length(acf),1));
      set(h,'color',[1 0 0]);
      set(h1,'color',[1 0 0]);
%%


subplot(2,2,3)
histogram(diff(swr{day}{epoch}{1}.starttime), 100, 'Normalization', 'cdf',...
    'FaceColor', [0 0 .3], 'LineStyle', ':');
title('Inter-SWR Interval CDF')
xlabel('Log10 time s')
ylabel('Probability')
set(gca, 'xscale', 'log')
set(gca, 'xtick', [0 2 4 6 8 10 20 40])
set(gca, 'ytick', 0:.2:1)
% axis tight
xlim([0 60])
ylim([0 1])