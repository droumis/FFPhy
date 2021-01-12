
%{
using position, dio, and manual review, determine rewarded visit intervals
- start: arrival at well = first lick (1 second later dio out if rewarded)
- end: animal leaves well area

%}
pconf = paramconfig;
getIngestionIntervals = 1;

plot_IngestionIntervals_pDay = 1;
show_figs = 1;
pause_figs = 0;
save_figs = 1;
save_fig_as = {'pdf', 'png'};

animals = {'JZ4', 'D10'};
day = 1;
epochs = [2,4];

%%
if getIngestionIntervals
    F = struct;
    for ian = 1:length(animals)
        animal = animals{ian};
        % load position
        p = load_data('filterframework', 'pos01', animal,'animpos', 0);
        posfields = strsplit(p.pos{day}{epoch}.fields);
        posidx = getFieldIndex(posfields, ...
            {'time' 'x-loess' 'y-loess' 'vel-loess'});
        
        % load dio
        dio = load_data('filterframework','DIO01',animal,'animpos',0);
        dio = dio.DIO;
        task = load_data('filterframework','task01',animal,'animpos',0);
        task = task.task;
        for e = 1:length(epochs)
            pos = p.pos{day}{epochs(e)}.data(:, posidx);
            % get times of reward dio output
            rewT = get_rewardTimes(task, dio, day, epochs(e));
            arrivs = [rewT(:,1)-1.100 rewT(:,2)];
            
            % get X,Y position of reward wells
            %             linpos = load_data('filterframework','linpos01',animal,'animpos',0);
            %             linpos = linpos.linpos;
            %             wellCoords = linpos{day}{epochs(e)}.wellSegmentInfo.wellCoord;
            
            % get times when animal departs well area, and when starting to move
            %             departs = [];
            %             posXYdepart = [];
            leaving = [];
            %             wellCoords = sort(wellCoords);
            %             departDistThresh = 30;
            leavingVelThresh = 4;
            for i = 1:size(rewT, 1)
                posIdxAtRewT = lookup(rewT(i,1), pos(:,1));
                velAtWellPosIdx = find(pos(posIdxAtRewT:end,4) > ...
                    leavingVelThresh, 1);
                movingAwayPosIdx = velAtWellPosIdx + posIdxAtRewT;
                try
                    leaving = [leaving; pos(movingAwayPosIdx, 1)];
                catch
                    fprintf('out of bounds\n');
                    arrivs(i,:) = [];
                    continue
                end
                %             % use x coord of position to determine which well he is at
                %             atWell = knnsearch(wellCoords(:,1), pos(posIdxAtRewT,2));
                %             atWellCoord = wellCoords(atWell, :);
                %             yposFarFromWell = find(pos(posIdxAtRewT:end,3) - atWellCoord(2) ...
                %                 > departDistThresh,1);
                %             departPosIdx = yposFarFromWell + posIdxAtRewT;
                %             departs = [departs; pos(departPosIdx, 1)];
                %         posXYdepart = [posXYdepart; pos(departPosIdx, 1)];
            end
            
            F(ian).ingIntervals{e} = [arrivs(:,1) leaving];
            F(ian).pos{e} = pos;
        end
        F(ian).animal = animal;
    end
    % get lick bout intervals
    % [~, boutTimes] = getLickBoutLicks(animal, [repmat(day, length(eps),1) eps']);
    % boutTimesall = vertcat(boutTimes{day}{eps});
end
%% plot pos, rew intervals, rew dio output,
if plot_IngestionIntervals_pDay
    label = 'ingestionIntervals';
    figname = sprintf('%s-pDay', label);
    Pp=load_plotting_params({'defaults', figname});
    for ian = 1:length(animals)
        animal = F(ian).animal;
        ifig = init_plot(show_figs, Pp.position);
        for e = 1:length(epochs)
            % collect epoch data
            pos = F(ian).pos{e};
            ingIntervals = F(ian).ingIntervals{e};
            
            subaxis(2,1,e,Pp.posparams{:});
            plot(pos(:,1), pos(:,[2 3 4]))
            axis tight
            ylim([0 220])
            xlabel('time')
            
            % plot ingestion intervals
            line(repmat(ingIntervals(:,1),1,2)', ...
                repmat(ylim, length(ingIntervals(:,1)), 1)', ...
                'linewidth', 2, 'color', 'k')
            
            line(repmat(ingIntervals(:,2),1,2)', ...
                repmat(ylim, length(ingIntervals(:,2)), 1)', ...
                'linewidth', 2, 'color', 'r')
        end
        % super axis
        stit = sprintf('%s %s', figname, animal);
        setSuperAxTitle(stit);
        if pause_figs
            pause
        end
        if save_figs
            strsave = save_figure([pconf.andef{4} '/' figname],...
                stit, 'savefigas', save_fig_as);
        end
    end
end
%% plot pos X,Y and well coordinates
% cla
% figure(2)
% plot(pos(:,2), pos(:,3))
% scatter(wellCoords(:,1),wellCoords(:,2), 5000, 'r.')
























