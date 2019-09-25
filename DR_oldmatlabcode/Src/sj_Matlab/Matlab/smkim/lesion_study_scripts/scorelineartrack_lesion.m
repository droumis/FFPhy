
%{
clear('summary_lineartrack');

% count various stuffs for learning curves
summary_lineartrack.labels = { 'M02' 'M03' 'M06' 'M12' 'M13' 'M14' 'M16' 'M17' 'M19' 'M20' 'M22' 'M24' 'M25' 'M26' };
summary_lineartrack.islesion = [  1     1     0     1     1     1     1     1     1     1     1     0     0     0    ];
summary_lineartrack.days = 1:2;

% summary_lineartrack variables
zeroblock = zeros([numel(summary_lineartrack.days) numel(summary_lineartrack.labels)]);
summary_lineartrack.total = zeroblock;
summary_lineartrack.correct = zeroblock;
summary_lineartrack.sinuosity = zeroblock;
summary_lineartrack.runningspeed = zeroblock;
summary_lineartrack.immobilityfraction = zeroblock;
summary_lineartrack.pathlength = zeroblock;
summary_lineartrack.tripduration = zeroblock;

for s = 1:numel(summary_lineartrack.labels)
    load([summary_lineartrack.labels{s} '_lineartrack_trajectories.mat']);
    for d = summary_lineartrack.days

        % collapse sessions 1 and 2 within each day
        endzone = [ trajectories{d}{1}.endzone trajectories{d}{2}.endzone ];
        startzone = [ trajectories{d}{1}.startzone trajectories{d}{2}.startzone ];
        correct = [ trajectories{d}{1}.correct trajectories{d}{2}.correct ];
        smoothedpos_snippets = [ trajectories{d}{1}.smoothedpos_snippets ...
          trajectories{d}{2}.smoothedpos_snippets ];

        % allocate temporary data containers
        runningspeed = [];
        pathlength = [];
        tripduration = [];
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
                sinuosity(end+1) = pathsinuosity(smoothedpos_snippets(i),15);
                % path length
                pathlength(end+1) = sum(hypot( ...
                  diff(smoothedpos_snippets(i).data(:,2)),diff(smoothedpos_snippets(i).data(:,3))));
                % trip duration
                tripduration(end+1) = ...
                  smoothedpos_snippets(i).data(end,1) - smoothedpos_snippets(i).data(1,1);
                % total count
                summary_lineartrack.total(d,s) = summary_lineartrack.total(d,s) + 1;
                summary_lineartrack.correct(d,s) = summary_lineartrack.correct(d,s) + correct(i);
            end
        end
        summary_lineartrack.sinuosity(d,s) = mean(sinuosity);
        summary_lineartrack.runningspeed(d,s) = mean(runningspeed);
        summary_lineartrack.immobilityfraction(d,s) = mean(immobilityfraction);
        summary_lineartrack.pathlength(d,s) = mean(pathlength);
        summary_lineartrack.tripduration(d,s) = mean(tripduration);
    end
    clear('trajectories');
end
save('/data15/smkim/lesion_study/behavior/summary_lineartrackperformance.mat','summary_lineartrack');
%}

d = summary_lineartrack.pathlength;
labelstring = 'path length between food wells (cm)';
yrange = [ ];
%
%%%%% mean and SEM plots
%{
d = summary_lineartrack.total;
labelstring = 'number of food well visits';
yrange = [ ];
%
d = summary_lineartrack.correct ./ summary_lineartrack.total;
labelstring = 'proportion correct visits';
yrange = [0.5 1.3];
%
d = summary_lineartrack.runningspeed;
labelstring = 'average running speed (cm/s)';
yrange = [ 0 80];
%
d = summary_lineartrack.immobilityfraction;
labelstring = 'fraction of time spent immobile';
yrange = [0 0.1];
%
d = summary_lineartrack.sinuosity;
labelstring = 'sinuosity of running path';
yrange = [0 0.2];
%
d = summary_lineartrack.tripduration;
labelstring = 'journey time between food wells (s)';
yrange = [0 10];
%
d = summary_lineartrack.pathlength;
labelstring = 'path length between food wells (cm)';
yrange = [ ];
%}

y_sham = mean(d(:,summary_lineartrack.islesion==0),2);
e_sham = std(d(:,summary_lineartrack.islesion==0),0,2)/sqrt(nnz(~summary_lineartrack.islesion));
y_lesion = mean(d(:,summary_lineartrack.islesion==1),2);
e_lesion = std(d(:,summary_lineartrack.islesion==1),0,2)/sqrt(nnz(summary_lineartrack.islesion));

figure(1);
set(gca,'XLim',[0.5 2.5],'FontSize',20,'XTick',summary_lineartrack.days,'LineWidth',5, ...
    'TickDir','out','TickLength',[0.03 0.03],'XTickLabel',{'pre','post'});
if ~isempty(yrange)
    set(gca,'YLim',yrange);
end
set(get(gca,'XLabel'),'FontSize',20,'Color','k','String','day');
set(get(gca,'YLabel'),'FontSize',20,'Color','k','String',labelstring);
hold on;
errorbar(summary_lineartrack.days,y_sham,e_sham, ...
  'Color','g','LineStyle','-','Marker','s','MarkerFaceColor','g', ...
  'LineWidth',4,'MarkerSize',16);
errorbar(summary_lineartrack.days,y_lesion,e_lesion, ...
  'Color','m','LineStyle','--','Marker','o','MarkerFaceColor','w', ...
  'LineWidth',4,'MarkerSize',16);
legend_handle = legend(['sham-operated control (' num2str(nnz(~summary_lineartrack.islesion)) ')'], ...
  ['lesion (' num2str(nnz(summary_lineartrack.islesion)) ')'], ...
    'Location','Best');
set(legend_handle,'Box','off','TextColor','k');

%{
%%%%% plot all
figure(2)
line(repmat((summary_lineartrack.days)',[1 nnz(summary_lineartrack.islesion==0)]), ...
  d(:,summary_lineartrack.islesion==0), ...
  'Color','g','LineStyle','-','Marker','s','MarkerFaceColor','w', ...
  'LineWidth',4,'MarkerSize',16);
line(repmat((summary_lineartrack.days)',[1 nnz(summary_lineartrack.islesion==1)]), ...
  d(:,summary_lineartrack.islesion==1), ...
  'Color','m','LineStyle','--','Marker','x','MarkerFaceColor','k', ...
  'LineWidth',4,'MarkerSize',16);
%
set(gca,'XLim',[0.5 10.5],'FontSize',20,'XTick',summary_lineartrack.days,'LineWidth',5, ...
    'TickDir','out','TickLength',[0.03 0.03]);
if ~isempty(yrange)
    set(gca,'YLim',yrange);
end
set(get(gca,'XLabel'),'FontSize',24,'Color','k','String','day');
set(get(gca,'YLabel'),'FontSize',20,'Color','k','String',labelstring);

%}

