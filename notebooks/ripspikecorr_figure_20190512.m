
%
loaddata = 0;
plotfigs = 1;
plot_heatmaps = 0;
pausefigs = 0;
savefigs = 1;
env = 'wtrack';

animals = {'JZ3'};
me = animaldef('demetris');

%% Loaddata add the r and p to the animal level
if loaddata
    paths = make_paths('paircorrVperformance', env);
    data = [];
    for an = 1:length(animals)
        data{an} = load_filter_output(paths.filtOutputDirectory, paths.filenamesave, ...
            'animal', animals{an});
        data{an}.results = [[data{an}.output.ntA]' [data{an}.output.ntB]' ...
            [data{an}.output.xc_perf_r]' [data{an}.output.xc_perf_p]' ...
            [data{an}.output.xc_diff_r]' [data{an}.output.xc_diff_p]'];
        data{an}.results_fields = ['ntA ntB perf_r perf_p diff_r diff_p'];
    end
end

%% Plot full figs

if plotfigs
    numrows = 8;
    binsize = .002;
    smoothing_length = 10;   % std of gaussian (in ms) used to smooth rip psth
    smoothing_width = round(smoothing_length*.001/binsize);   % smoothing width in number of bins
    kernel = gaussian(smoothing_width,smoothing_width*8);
    Pp = load_plotting_params('paircorrVperformance');
    paths = make_paths('paircorrVperformance', env);
    for an = 1:length(animals)
        animal = data{an}.animal{3};
        pairs = [[data{an}.output.ntA]' [data{an}.output.ntB]'];
        pz = cell2mat({data{an}.output.xc_diff_p}');
        pzsig = find(pz<Pp.alpha);
        ppairs = pairs(pzsig,:);
        for p = 1:length(ppairs(:,1))
            if savefigs && ~pausefigs;
                close all
                ifig = figure('Visible','off','units','normalized','position', ...
                    Pp.position);
            else
                ifig = figure('units','normalized','position',Pp.position);
            end
            set(gcf,'color','white')
            
            ntA = ppairs(p,1);
            ntB = ppairs(p,2);
            ppdataind = find(ismember(pairs, [ntA ntB],'rows'));
            idata = data{an}.output(ppdataind);
            ipair_ids = cell2mat({idata.pair_data.index}');
            days = unique(ipair_ids(:,1),'stable');
            numdays = length(days);
            
            %% plot ntA riptrig spiking (rows 1-2)
            for nt = 1:2
                rows = [1 3];
                for d = 1:numdays
                    day = days(d);
                    day_inds = find(ipair_ids(:,1)==day);
                    % combine the pairdata across epochs for each day and
                    % flip to make time go up
                    if nt == 1
                        dayspikes = flipud(cell2mat({idata.pair_data(day_inds).ntA_spike_raster}'));
                    else
                        dayspikes = flipud(cell2mat({idata.pair_data(day_inds).ntB_spike_raster}'));
                    end
                    [sx, sy] = find(dayspikes');
                    subplot(numrows,numdays,d+numdays*(rows(nt)-1))
                    f1 = scatter(sx,sy,3, 'filled');
%                     set(gca,'YDir','reverse')
                    f1.MarkerFaceAlpha = 0.15;
                    f1.MarkerFaceColor = [0 0 0];
                    axis tight
                    set(gca,'XTick',[]);
                    set(gca,'YTick',[]);
                    if d == 1;
                        ylabel({[sprintf('nt %d',ppairs(p,nt))];['swr #']})
                    end
                    subplot(numrows,numdays,d+numdays*rows(nt))
%                     t =idata.pair_data(day_inds(1)).time;
                        % Create smoothed PSTH of mean firing rate
                    numrips = length(dayspikes(:,1));
                    binsize = .002; %2ms currently
                    smoothing_length = 10;   % std of gaussian (in ms) used to smooth rip psth
                    smoothing_width = round(smoothing_length*.001/binsize);   % smoothing width in number of bins
                    kernel = gaussian(smoothing_width,smoothing_width*8);
                    smoothedpsth = smoothvect(sum(full(dayspikes),1)./(binsize*numrips),kernel);
                    plmin=min(smoothedpsth(5:end-5));
                    lh = area(smoothedpsth, plmin);
                    lh.LineWidth = 0.0001;
                    lh.FaceColor = [0 0 0];
                    try
                        ylim([plmin max(smoothedpsth)])
                        xlim([1 length(smoothedpsth)])
                    catch
                        fprintf('no spikes\n')
                        continue
                    end
                    set(gca,'XTick',[]);
                    set(gca,'YTick',[]);
%                     clear ylim;

%                     axis tight
                    % fuck ok i forgot that the rasters are binned into 2
                    % ms bins bc that how i did the xcorr.. so 5 bins is 10
                    % ms
%                     histogram(sx, 100, 'DisplayStyle', 'stairs','Normalization','probability'); %/numrips/binsize*1000;
%                     area([1 2 3],'k')
%                     axis tight
                    % the bar funciton is broken
                    
                    %                     % ripple line
                    %                     lx = round(max(sx)/2);
                    %                     ly = size(di.psth,1);
                    %                     line([lx lx],[1 ly],'Color',[0 0 0], 'LineStyle', '--')
                    %                     set(gca, 'XTickLabel',[])
                    %                     set(gca, 'YTickLabel',[])
                end
            end
                                

%             nbins=20;
%             histogram(sx,nbins, 'FaceColor', [.8 .8 .8], 'FaceAlpha', .1,...
%                 'DisplayStyle', 'bar', 'EdgeAlpha', .3, 'Normalization','probability');
%             Create smoothed PSTH of mean firing rate
%             smoothedpsth = smoothvect(sum(dayspikes),1)./(binsize*data.noevents),kernel);
%             wp line assumes there are 2 wtrack epochs per day
%                                 firstepochlen = length(idata.pair_data(day_inds(1)).ntA_spike_raster);
%                                 line([1 500],[firstepochlen firstepochlen])
            %
            %                     subplot(7,numdays,d+numdays)
                                
            % % %                     kern = gaussian();
            %                                 line([lx lx],[1 max(h.Values)],'Color',[0 0 0], 'LineStyle', '--')
            %                     zhist = smoothvec(zscore(sum(xi)),kern);
            %                     line(zhist, 'filled')
            %                 end
            %                 % plot aggregate for each day
            
            %% plot ntAB pair corr per day (row 3)
            %           ec = [xc d ep riptime]
            try
                exc = cell2mat({idata.pair_data.ntAB_excesscorr}');
            catch
                exc = [idata.pair_data.ntAB_excesscorr]';
            end
            riptimes = cell2mat({idata.pair_data.ripstarttimes}');
            
            daylen = cellfun(@(x) length(x),{idata.pair_data.ntAB_excesscorr}',...
                'un', 1);
            dayids = cellfun(@(x) x(1),{idata.pair_data.index}', 'un', 1);
            daymat = cell2mat(arrayfun(@(x,y) repmat(x,y,1),dayids,daylen, ...
                'un',0));
            
            try
                perf = cell2mat({idata.pair_data.performance}');
                perfdiff = cell2mat({idata.pair_data.performance_diff}');
            catch
                perf = [idata.pair_data.performance]';
                perfdiff = [idata.pair_data.performance_diff]';
            end
            
            result_fields = ['excesscorr ripstarttime day perf perfdiff'];
            ipair_result_mat = [exc riptimes daymat perf perfdiff];
            
%             subplot(numrows,numdays,[numdays*4+1:numdays*5])
%             boxplot(ipair_result_mat(:,1),ipair_result_mat(:,3),'Notch','on')
%             ylabel('spike corr')
%             xlabel('day')
            % instead, plot over time, for each day
            for d = 1:numdays
                day = days(d);
                subplot(numrows,numdays,numdays*4+d)
                y = ipair_result_mat(ipair_result_mat(:,3)==day, 1);
%                 x = ipair_result_mat(ipair_result_mat(:,3)==day, 2);
                plot(y,'.b')
                ylim([min(ipair_result_mat(:,1)),max(ipair_result_mat(:,1))])
                xlim([1 length(y)])
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
            end
            
            %% plot performance and diff (row 4)
            for d = 1:numdays
                day = days(d);
                subplot(numrows,numdays,numdays*5+d)
                y = ipair_result_mat(ipair_result_mat(:,3)==day, 4);
                y = y(~isnan(y));
                plot(y,'b')
                ylim([min(ipair_result_mat(:,4)),max(ipair_result_mat(:,4))])
                xlim([1 length(y)])
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
                subplot(numrows,numdays,numdays*6+d)
                dy = ipair_result_mat(ipair_result_mat(:,3)==day, 5);
                dy = dy(~isnan(dy));
                plot(dy,'r')
                ylim([min(ipair_result_mat(:,5)),max(ipair_result_mat(:,5))])
                xlim([1 length(dy)])
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
            end
%             subplot(numrows,numdays,[numdays*5+1:numdays*6])
%             plot(cell2mat({idata.pair_data.performance}'), '.')
%             hold on
%             plot(cell2mat({idata.pair_data.performance_diff}')*20+.5,'.')
%             line([cumsum(daylen)'; cumsum(daylen)'], [0 1], 'Color', [.9 .9 .9], ...
%                 'LineStyle','-', 'LineWidth', 1)
%             axis tight
%             hold off
%             xlabel('trials')
%             ylabel('perf+diff')
            
            %             modecol = find(cell2mat(cellfun(@(x) strcmp(x,'mode'), strsplit(statespace.pcFields, ' '), 'UniformOutput', false)));
            %             lower5col = find(cell2mat(cellfun(@(x) strcmp(x,'lower5'), strsplit(statespace.pcFields, ' '), 'UniformOutput', false)));
            %             upper5col = find(cell2mat(cellfun(@(x) strcmp(x,'upper5'), strsplit(statespace.pcFields, ' '), 'UniformOutput', false)));
            %             certaintycol = find(cell2mat(cellfun(@(x) strcmp(x,'certainty'), strsplit(statespace.pcFields, ' '), 'UniformOutput', false)));
            %
            %             line([cumsum(statespace.eplengths); cumsum(statespace.eplengths)], [0 1], 'Color', [.9 .9 .9], 'LineStyle','-', 'LineWidth', 1)
            %             plot(tALL, pc(2:end,modecol),'b-', 'LineWidth', 2); %plot behavior SS score
            %             errfillAll = fill([tALL fliplr(tALL)],[pc(2:end,lower5col); flipud(pc(2:end,upper5col))],[0 0 1],'linestyle','none');
            %             set(errfillAll, 'FaceAlpha', .2)
            %             hold on; [x, y] = find(behavperform > 0);
            %             h = plot(x,y-0.03,'+'); set(h, 'MarkerFaceColor','none');
            %             set(h, 'MarkerEdgeColor', [.2 .7 .2]);
            %             hold on; [x, y] = find(behavperform == 0);
            %             h = plot(x,y-0.05,'+'); set(h, 'MarkerFaceColor', 'none');
            %             set(h, 'MarkerEdgeColor', [.5 .5 .5]);
            %             axis([1 tALL(end)  0 1.05]);
            %             line([1 tALL(end)], [chance  chance], 'Color', [.5 .5 .5], 'LineStyle', '--');
            %             xlabel('Trial Number')
            %             ylabel([{'Probability of a'};{'Correct Response'}])
            
            
            %% plot scatters (row 5)
            halfdays = round(numdays/2);
%             subplot(numrows,numdays,[numdays*7+1:(numdays*8-halfdays)])
            subplot(numrows,numdays,numdays*7+1)
            scatter(ipair_result_mat(:,1),ipair_result_mat(:,4), 4,ipair_result_mat(:,3))
            ylabel(sprintf('perf r:%.3f p:%.3f',idata.xc_perf_r, idata.xc_perf_p))
            xlabel('xc')
            axis tight
            lsline
%             subplot(5,numdays,[(numdays*8-halfdays+1):numdays*8])
            subplot(numrows,numdays,numdays*7+2)
            scatter(ipair_result_mat(:,1),ipair_result_mat(:,5), 4,ipair_result_mat(:,3))
            ylabel(sprintf('diff r:%.3f p:%.4f',idata.xc_diff_r, idata.xc_diff_p))
            xlabel('xc')
            axis tight
            lsline
            
            %% ---- super title and colorbar----
            sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
            sprtit = sprintf('riptrigmu xcVbs %s %s - nt%d nt%d', env, animal, ntA, ntB);
            iStitle = text(.5, .95, {sprtit}, 'Parent', sprtitleax, 'Units', ...
                'normalized');
            set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
                'horizontalAlignment', 'center');
            %% ---- pause, save figs ----
            if pausefigs
                pause
            end
            if savefigs
                save_figure(paths.figdirectory, paths.filenamesave, sprtit)
            end
        end
    end
end



%%
if plot_heatmaps
    paths = make_paths('paircorrVperformance_heatmap', env);
    Pp = load_plotting_params('paircorrVperformance_hm');
    for an = 1:length(animals)
        if savefigs && ~pausefigs;
            close all
            ifig = figure('Visible','off','units','normalized','position', ...
                Pp.position);
        else
            %             clf
            ifig = figure('units','normalized','position',Pp.position);
        end
        set(gcf,'color','white')
        animal = data{an}.animal{3};
        nta = unique(data{an}.results(:,1),'stable');
        ntb = unique(data{an}.results(:,2),'stable');
        %         o = nan(length(xs),length(ys));
        %         for i = 1:length(data{an}.results(:,1))
        %             w = data{an}.results(i,:);
        %             x = round(w(1));
        %             y = round(w(2));
        %             o(x,y) = w(6);
        %         end
        Data2D = accumarray([data{an}.results(:,1), data{an}.results(:,2)], ...
            data{an}.results(:,6),[max(nta) max(ntb)],[],nan);
        s = sort(Data2D(:), 'ascend');
        heatmap_custom(Data2D, 1:max(ntb), 1:max(nta), '%0.3f', 'TickAngle', 0,...
            'ShowAllTicks', true,'Colormap', 'cool', 'Colorbar', true, ...
            'UseLogColormap', false, 'NaNColor', [0 0 0], 'MinColorValue', s(1),...
            'MaxColorValue', s(40), 'GridLines', ':', 'FontSize', 4, 'TextColor', 'k');
        ylabel('nta')
        xlabel('ntb')
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontName','arial','fontsize',18)
        %% ---- super title and colorbar----
        sprtitleax = axes('Position',[0 0 1 1],'Visible','off', 'Parent', ifig);
        sprtit = sprintf('%s MEC %s ripcorrVperfDIFF maxP e%.0f',animal,env,log(s(40)));
        iStitle = text(.5, .95, {sprtit}, 'Parent', sprtitleax, 'Units', ...
            'normalized');
        set(iStitle,'FontWeight','bold','Color','k', 'FontName', 'Arial', ...
            'horizontalAlignment', 'center');
        %% ---- pause, save figs ----
        %% ---- pause, save figs ----
        if pausefigs
            pause
        end
        if savefigs
            save_figure(paths.figdirectory, paths.filenamesave, sprtit)
        end
    end
end




































































































































