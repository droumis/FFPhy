
%{
clear('summary_Wtrack');

% count various stuffs for learning curves
summary_Wtrack.labels = { 'M02' 'M03' 'M06' 'M12' 'M13' 'M14' 'M16' 'M17' 'M19' 'M20' 'M22' 'M24' 'M25' 'M26' };
summary_Wtrack.islesion = [  1     1     0     1     1     1     1     1     1     1     1     0     0     0    ];
summary_Wtrack.days = 1:10;

% summary_Wtrack variables
zeroblock = zeros([numel(summary_Wtrack.days) numel(summary_Wtrack.labels)]);
summary_Wtrack.total = zeroblock;
summary_Wtrack.correct = zeroblock;
summary_Wtrack.total_fromside = zeroblock;
summary_Wtrack.correct_fromside = zeroblock;
summary_Wtrack.total_fromcenter = zeroblock;
summary_Wtrack.correct_fromcenter = zeroblock;
summary_Wtrack.error_revisit = zeroblock;
summary_Wtrack.right = zeroblock;
summary_Wtrack.center = zeroblock;
summary_Wtrack.left = zeroblock;
summary_Wtrack.sinuosity = zeroblock;
summary_Wtrack.runningspeed = zeroblock;
summary_Wtrack.immobilityfraction = zeroblock;
summary_Wtrack.pathlength_side_side = zeroblock;
summary_Wtrack.pathlength_side_center = zeroblock;
summary_Wtrack.tripduration_side_side = zeroblock;
summary_Wtrack.tripduration_side_center = zeroblock;

for s = 1:numel(summary_Wtrack.labels)
    load([summary_Wtrack.labels{s} '_Wtrack_trajectories.mat']);
    for d = summary_Wtrack.days

        % collapse sessions 1 and 2 within each day
        endzone = [ trajectories{d}{1}.endzone trajectories{d}{2}.endzone ];
        startzone = [ trajectories{d}{1}.startzone trajectories{d}{2}.startzone ];
        correct = [ trajectories{d}{1}.correct trajectories{d}{2}.correct ];
        smoothedpos_snippets = [ trajectories{d}{1}.smoothedpos_snippets ...
          trajectories{d}{2}.smoothedpos_snippets ];

        % allocate temporary data containers
        runningspeed = [];
        pathlength_side_side = [];
        pathlength_side_center = [];
        tripduration_side_side = [];
        tripduration_side_center = [];
        immobilityfraction = [];
        sinuosity = [];

        for i = 1:length(startzone);
            % consider only traj with well-defined startzone and endzone
            if ~isnan(startzone(i)) & ~isnan(endzone(i))
                % running speed
                runningspeed(end+1) = mean(hypot(smoothedpos_snippets(i).data(:,4), ...
                  smoothedpos_snippets(i).data(:,5)));
                % immobility fraction
                immobilityfraction(end+1) = nnz(smoothedpos_snippets(i).data(:,6))/ ...
                  size(smoothedpos_snippets(i).data,1);
                % sinuosity
                sinuosity(end+1) = pathsinuosity(smoothedpos_snippets(i),10);
                % total count
                summary_Wtrack.total(d,s) = summary_Wtrack.total(d,s) + 1;
                summary_Wtrack.correct(d,s) = summary_Wtrack.correct(d,s) + correct(i);
                % score correct/error on inbound versus outbound
                if (startzone(i) == 1) | (startzone(i) == 3)
                    summary_Wtrack.total_fromside(d,s) = summary_Wtrack.total_fromside(d,s) + 1;
                    if endzone(i) == 2
                        pathlength_side_center(end+1) = sum(hypot( ...
                          diff(smoothedpos_snippets(i).data(:,2)),diff(smoothedpos_snippets(i).data(:,3))));
                        tripduration_side_center(end+1) = ...
                          smoothedpos_snippets(i).data(end,1) - smoothedpos_snippets(i).data(1,1);
                    elseif endzone(i) ~= startzone(i)
                        pathlength_side_side(end+1) = sum(hypot( ...
                          diff(smoothedpos_snippets(i).data(:,2)),diff(smoothedpos_snippets(i).data(:,3))));
                        tripduration_side_side(end+1) = ...
                          smoothedpos_snippets(i).data(end,1) - smoothedpos_snippets(i).data(1,1);
                    end
                    if correct(i)
                        summary_Wtrack.correct_fromside(d,s) = summary_Wtrack.correct_fromside(d,s) + 1;
                    end
                elseif startzone(i) == 2
                    summary_Wtrack.total_fromcenter(d,s) = summary_Wtrack.total_fromcenter(d,s) + 1;
                    if endzone(i) ~= 2
                        pathlength_side_center(end+1) = sum(hypot( ...
                          diff(smoothedpos_snippets(i).data(:,2)),diff(smoothedpos_snippets(i).data(:,3))));
                        tripduration_side_center(end+1) = ...
                          smoothedpos_snippets(i).data(end,1) - smoothedpos_snippets(i).data(1,1);
                    end
                    if correct(i)
                        summary_Wtrack.correct_fromcenter(d,s) = summary_Wtrack.correct_fromcenter(d,s) + 1;
                    end
                end
                % score perseverative error
                if startzone(i) == endzone(i)
                    summary_Wtrack.error_revisit(d,s) = summary_Wtrack.error_revisit(d,s) + 1;
                end
                % score left versus right
                if endzone(i) == 1
                    summary_Wtrack.right(d,s) = summary_Wtrack.right(d,s) + 1;
                elseif endzone(i) == 2
                    summary_Wtrack.center(d,s) = summary_Wtrack.center(d,s) + 1;
                elseif endzone(i) == 3
                    summary_Wtrack.left(d,s) = summary_Wtrack.left(d,s) + 1;
                end 
            end
        end
        summary_Wtrack.sinuosity(d,s) = mean(sinuosity);
        summary_Wtrack.runningspeed(d,s) = mean(runningspeed);
        summary_Wtrack.immobilityfraction(d,s) = mean(immobilityfraction);
        summary_Wtrack.pathlength_side_side(d,s) = mean(pathlength_side_side);
        summary_Wtrack.pathlength_side_center(d,s) = mean(pathlength_side_center);
        summary_Wtrack.tripduration_side_side(d,s) = mean(tripduration_side_side);
        summary_Wtrack.tripduration_side_center(d,s) = mean(tripduration_side_center);
    end
    clear('trajectories');
end
save('/data15/smkim/lesion_study/behavior/summary_Wtrackperformance.mat','summary_Wtrack');
%}

%
%%%%% mean and SEM plots


d = summary_Wtrack.sinuosity;
labelstring = 'sinuosity of running path';
yrange = [0 0.4];
%{
d = summary_Wtrack.total;
labelstring = 'number of food well visits';
yrange = [ ];
%
d = summary_Wtrack.correct ./ summary_Wtrack.total;
labelstring = 'proportion correct visits';
yrange = [0 1];
%
d = summary_Wtrack.correct_fromcenter ./ summary_Wtrack.total_fromcenter;
labelstring = 'correct choices from center arm';
yrange = [0 1];
%
d = summary_Wtrack.correct_fromside ./ summary_Wtrack.total_fromside;
labelstring = 'correct choices from side arms';
yrange = [0 1];
%
d = summary_Wtrack.error_revisit ./ summary_Wtrack.total;
labelstring = 'proportion perseverative errors';
yrange = [0 0.4];
%
d = (max(summary_Wtrack.right,summary_Wtrack.left) ./ (summary_Wtrack.right + summary_Wtrack.left) - 0.5);
yrange = [0 0.2 ];
labelstring = 'side bias'
%
d = summary_Wtrack.runningspeed;
labelstring = 'average running speed (cm/s)';
yrange = [ 0 80];
%
d = summary_Wtrack.immobilityfraction;
labelstring = 'fraction of time spent immobile';
yrange = [0 0.2];
%
d = summary_Wtrack.sinuosity;
labelstring = 'sinuosity of running path';
yrange = [0 0.4];
%
d = summary_Wtrack.tripduration_side_center;
labelstring = 'journey time between food wells (s)';
yrange = [ ];
%
d = summary_Wtrack.pathlength_side_center;
labelstring = 'path length between food wells (cm)';
yrange = [ ];
%}

y_sham = mean(d(:,summary_Wtrack.islesion==0),2);
e_sham = std(d(:,summary_Wtrack.islesion==0),0,2)/sqrt(nnz(~summary_Wtrack.islesion));
y_lesion = mean(d(:,summary_Wtrack.islesion==1),2);
e_lesion = std(d(:,summary_Wtrack.islesion==1),0,2)/sqrt(nnz(summary_Wtrack.islesion));

figure(1);
set(gca,'XLim',[0.5 10.5],'FontSize',20,'XTick',summary_Wtrack.days,'LineWidth',5, ...
    'TickDir','out','TickLength',[0.03 0.03]);
if ~isempty(yrange)
    set(gca,'YLim',yrange);
end
set(get(gca,'XLabel'),'FontSize',20,'Color','k','String','day');
set(get(gca,'YLabel'),'FontSize',20,'Color','k','String',labelstring);
hold on;
errorbar(summary_Wtrack.days,y_sham,e_sham, ...
  'Color','g','LineStyle','-','Marker','s','MarkerFaceColor','g', ...
  'LineWidth',4,'MarkerSize',16);
errorbar(summary_Wtrack.days,y_lesion,e_lesion, ...
  'Color','m','LineStyle','--','Marker','o','MarkerFaceColor','w', ...
  'LineWidth',4,'MarkerSize',16);
legend_handle = legend(['sham-operated control (' num2str(nnz(~summary_Wtrack.islesion)) ')'], ...
  ['lesion (' num2str(nnz(summary_Wtrack.islesion)) ')'], ...
    'Location','Best');
set(legend_handle,'Box','off','TextColor','k');

%{
%%%%% plot all
figure(2)
line(repmat((summary_Wtrack.days)',[1 nnz(summary_Wtrack.islesion==0)]), ...
  d(:,summary_Wtrack.islesion==0), ...
  'Color','g','LineStyle','-','Marker','s','MarkerFaceColor','w', ...
  'LineWidth',4,'MarkerSize',16);
line(repmat((summary_Wtrack.days)',[1 nnz(summary_Wtrack.islesion==1)]), ...
  d(:,summary_Wtrack.islesion==1), ...
  'Color','m','LineStyle','--','Marker','x','MarkerFaceColor','k', ...
  'LineWidth',4,'MarkerSize',16);
%
set(gca,'XLim',[0.5 10.5],'FontSize',20,'XTick',summary_Wtrack.days,'LineWidth',5, ...
    'TickDir','out','TickLength',[0.03 0.03]);
if ~isempty(yrange)
    set(gca,'YLim',yrange);
end
set(get(gca,'XLabel'),'FontSize',24,'Color','k','String','day');
set(get(gca,'YLabel'),'FontSize',20,'Color','k','String',labelstring);

%}

