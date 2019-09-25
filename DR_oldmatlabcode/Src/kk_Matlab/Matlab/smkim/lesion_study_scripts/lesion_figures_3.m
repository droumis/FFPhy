

%{
% separate analysis of inbound and outbound trials on Wtrack
clear;
load('behavperform_2');
lesion_idx = find(strcmp({behavperform(:).group},'lesion'));
control_idx = find(strcmp({behavperform(:).group},'control'));

for s = 1:length(behavperform)
  if strcmp(behavperform(s).group,'control')
    subject_label = sprintf('Control subject #%d',find(s==control_idx));
  elseif strcmp(behavperform(s).group,'lesion')
    subject_label = sprintf('Hippocampal lesion subject #%d',find(s==lesion_idx));
  else
    error('i made a mistake');
  end
  
  % load variable 'journeys' into workspace
  load(sprintf('%s_Wtrack_journeys.mat',behavperform(s).subject));
  % classify journeys as inbound/outbound and collect information (keep only
  % those that have defined scores, e.g. not NaN)
  inbound_summary(s) = struct(...
      'subject',subject_label, ...
      'group',behavperform(s).group, ...
      'number_trials',nan(10,1), ...
      'proportion_correct',nan(10,1), ...
      'proportion_revisit_error',nan(10,1), ...
      'proportion_side_error',nan(10,1), ...
      'mean_running_speed',nan(10,1), ...
      'median_visit_duration',nan(10,1));
  outbound_summary(s) = struct(...
      'subject',subject_label, ...
      'group',behavperform(s).group, ...
      'number_trials',nan(10,1), ...
      'proportion_correct',nan(10,1), ...
      'proportion_revisit_error',nan(10,1), ...
      'proportion_side_error',nan(10,1), ...
      'mean_running_speed',nan(10,1), ...
      'median_visit_duration',nan(10,1));
  for d = 1:10
    inbound_trials = struct( ...
        'origin',{}, ...
        'destination',{}, ...
        'correct',{}, ...
        'visit_duration',{});
    outbound_trials = struct( ...
        'origin',{}, ...
        'destination',{}, ...
        'correct',{}, ...
        'visit_duration',{});
    inbound_speed = zeros(0,1);
    outbound_speed = zeros(0,1);
    for r = 1:2
      for j = 1:length(journeys{d}{r}.correct)
        if ~isnan(journeys{d}{r}.correct(j))
          if (journeys{d}{r}.startzone(j) == 2)
            outbound_trials(end+1) = struct( ...
              'origin',journeys{d}{r}.startzone(j), ...
              'destination',journeys{d}{r}.endzone(j), ...
              'correct',journeys{d}{r}.correct(j), ...
              'visit_duration',journeys{d}{r}.visitduration(j) );
            outbound_speed = [outbound_speed; hypot( ...
                journeys{d}{r}.smoothedpos_snippets(j).data(:,4), ...
                journeys{d}{r}.smoothedpos_snippets(j).data(:,5))];
          else
            inbound_trials(end+1) = struct( ...
              'origin',journeys{d}{r}.startzone(j), ...
              'destination',journeys{d}{r}.endzone(j), ...
              'correct',journeys{d}{r}.correct(j), ...
              'visit_duration',journeys{d}{r}.visitduration(j) );
            inbound_speed = [inbound_speed; hypot( ...
                journeys{d}{r}.smoothedpos_snippets(j).data(:,4), ...
                journeys{d}{r}.smoothedpos_snippets(j).data(:,5))];
          end
        else
          continue;
        end
      end
    end
    % inbound statistics
    inbound_summary(s).number_trials(d) = length(inbound_trials);
    inbound_summary(s).proportion_correct(d) = ...
        nnz([inbound_trials(:).correct]) / length(inbound_trials);
    inbound_summary(s).proportion_revisit_error(d) = ...
        nnz([inbound_trials(:).destination] == [inbound_trials(:).origin]) / ...
        length(inbound_trials);
    inbound_summary(s).proportion_side_error(d) = ...
        nnz(~[inbound_trials(:).correct] & ...
        ([inbound_trials(:).origin] ~= [inbound_trials(:).destination])) / ...
        length(inbound_trials);
    inbound_summary(s).mean_running_speed(d) = mean(inbound_speed);
    v_in = [inbound_trials(:).visit_duration]; v_in(isnan(v_in)) = [];
    inbound_summary(s).median_visit_duration(d) = median(v_in);

    % outbound statistics
    outbound_summary(s).number_trials(d) = length(outbound_trials);
    outbound_summary(s).proportion_correct(d) = ...
        nnz([outbound_trials(:).correct]) / length(outbound_trials);
    outbound_summary(s).proportion_revisit_error(d) = ...
        nnz([outbound_trials(:).destination] == [outbound_trials(:).origin]) / ...
        length(outbound_trials);
    outbound_summary(s).proportion_side_error(d) = ...
        nnz(~[outbound_trials(:).correct] & ...
        ([outbound_trials(:).origin] ~= [outbound_trials(:).destination])) / ...
        length(outbound_trials);
    outbound_summary(s).mean_running_speed(d) = mean(outbound_speed);
    v_out = [outbound_trials(:).visit_duration]; v_out(isnan(v_out)) = [];
    outbound_summary(s).median_visit_duration(d) = median(v_out);
  end
end

save('Wtrack_inbound_summary.mat','inbound_summary');
save('Wtrack_outbound_summary.mat','outbound_summary');
%}

%
% summary median+range for W track, separated by inbound and outbound
clear;
load('Wtrack_inbound_summary');
load('Wtrack_outbound_summary');
lesion_idx = find(strcmp({inbound_summary(:).group},'lesion'));
control_idx = find(strcmp({inbound_summary(:).group},'control'));

FONTSIZE = 14;
axes_style.Color = 'none';
axes_style.XLim = [0.5 10.5];
axes_style.XTick = 1:10;
axes_style.FontSize = FONTSIZE;
axes_style.FontName = 'Arial';
axes_style.TickDir = 'out';
axes_style.TickLength = [0.025 0.025];
axes_style.Box = 'off';
axes_style.Color = 'none';
axes_style.LineWidth = 2;

lesion_style.Color = 'r';
lesion_style.Marker = 'o';
lesion_style.MarkerSize = 6;
lesion_style.MarkerFaceColor = 'none';
lesion_style.LineStyle = 'none';
lesion_style.LineWidth = 2;
control_style.Color = 'k';
control_style.Marker = 's';
control_style.MarkerSize = 6;
control_style.MarkerFaceColor = 'k';
control_style.LineStyle = 'none';
control_style.LineWidth = 2;

xdata = (axes_style.XTick)';
lesion_xoffset = +0.13;
control_xoffset = -0.13;

%trial_type = 'Inbound';
trial_type = 'Outbound';

if strcmp(trial_type,'Inbound')
  correct = [inbound_summary(:).proportion_correct];
  total = [inbound_summary(:).number_trials];
  speed = [inbound_summary(:).mean_running_speed];
  duration = [inbound_summary(:).median_visit_duration];
elseif strcmp(trial_type,'Outbound');
  correct = [outbound_summary(:).proportion_correct];
  total = [outbound_summary(:).number_trials];
  speed = [outbound_summary(:).mean_running_speed];
  duration = [outbound_summary(:).median_visit_duration];
else
  error('trial_type must be Inbound or Outbound');
end

%ylim = [-0.0333 1];
%ax = axes(axes_style,'YLim',ylim,'YTick',0:0.25:1, ...
%    'Position',[0.15 0.5 0.55 0.35]);
%text(mean(axes_style.XLim),1.1*(ylim(2) - ylim(1)) + ylim(1), ...
%    trial_type,'VerticalAlignment','bottom', ...
%    'HorizontalAlignment','center','FontSize',FONTSIZE);
%set(get(ax,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
%    'String','Proportion correct trials');
%set(get(ax,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
%    'String','Day');
%% work-around for NaN values
%for d = 1:length(xdata)
%  tmp = correct(d,lesion_idx);
%  tmp(isnan(tmp)) = [];
%  z(d) = median(tmp);
%end
%line(xdata,z,'LineStyle','-','LineWidth',6, ...
%    'Color',lesion_style.Color);
%line(xdata,median(correct(:,control_idx),2),'LineStyle','-','LineWidth',6, ...
%    'Color',control_style.Color);
%line(xdata + lesion_xoffset,correct(:,lesion_idx),lesion_style);
%line(xdata + control_xoffset,correct(:,control_idx),control_style);
%pause;
%print(gcf,'-depsc','-r600', ...
%    sprintf('/home/smkim/lesion_manuscript/revisions/Wtrack_correct_%s.eps',trial_type));
%delete(gcf);

ylim = [-5 150];
ax = axes(axes_style,'YLim',ylim,'YTick',0:50:250, ...
    'Position',[0.15 0.5 0.55 0.3]);
text(mean(axes_style.XLim),1.1*(ylim(2) - ylim(1)) + ylim(1), ...
    trial_type,'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center','FontSize',FONTSIZE);
set(get(ax,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Number of trials performed');
set(get(ax,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Day');
% work-around for NaN values
for d = 1:length(xdata)
  tmp = total(d,lesion_idx);
  tmp(isnan(tmp)) = [];
  z(d) = median(tmp);
end
line(xdata,z,'LineStyle','-','LineWidth',6, ...
    'Color',lesion_style.Color);
line(xdata,median(total(:,control_idx),2),'LineStyle','-','LineWidth',6, ...
    'Color',control_style.Color);
line(xdata + lesion_xoffset,total(:,lesion_idx),lesion_style);
line(xdata + control_xoffset,total(:,control_idx),control_style);
pause;
print(gcf,'-depsc','-r600', ...
    sprintf('/home/smkim/lesion_manuscript/revisions/Wtrack_total_%s.eps',trial_type));
delete(gcf);

ylim = [0 100];
ax = axes(axes_style,'YLim',ylim,'YTick',0:20:100, ...
    'Position',[0.15 0.5 0.55 0.3]);
text(mean(axes_style.XLim),1.1*(ylim(2) - ylim(1)) + ylim(1), ...
    trial_type,'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center','FontSize',FONTSIZE);
set(get(ax,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Running speed (cm/s)');
set(get(ax,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Day');
% work-around for NaN values
for d = 1:length(xdata)
  tmp = speed(d,lesion_idx);
  tmp(isnan(tmp)) = [];
  z(d) = median(tmp);
end
line(xdata,z,'LineStyle','-','LineWidth',6, ...
    'Color',lesion_style.Color);
line(xdata,median(speed(:,control_idx),2),'LineStyle','-','LineWidth',6, ...
    'Color',control_style.Color);
line(xdata + lesion_xoffset,speed(:,lesion_idx),lesion_style);
line(xdata + control_xoffset,speed(:,control_idx),control_style);
pause;
print(gcf,'-depsc','-r600', ...
    sprintf('/home/smkim/lesion_manuscript/revisions/Wtrack_speed_%s.eps',trial_type));
delete(gcf);

ylim = [0 20];
ax = axes(axes_style,'YLim',ylim,'YTick',0:4:20, ...
    'Position',[0.15 0.5 0.55 0.3]);
text(mean(axes_style.XLim),1.1*(ylim(2) - ylim(1)) + ylim(1), ...
    trial_type,'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center','FontSize',FONTSIZE);
set(get(ax,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Average duration of food-well visit (s)');
set(get(ax,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Day');
% work-around for NaN values
for d = 1:length(xdata)
  tmp = duration(d,lesion_idx);
  tmp(isnan(tmp)) = [];
  z(d) = median(tmp);
end
line(xdata,z,'LineStyle','-','LineWidth',6, ...
    'Color',lesion_style.Color);
line(xdata,median(duration(:,control_idx),2),'LineStyle','-','LineWidth',6, ...
    'Color',control_style.Color);
line(xdata + lesion_xoffset,duration(:,lesion_idx),lesion_style);
line(xdata + control_xoffset,duration(:,control_idx),control_style);
pause;
print(gcf,'-depsc','-r600', ...
    sprintf('/home/smkim/lesion_manuscript/revisions/Wtrack_duration_%s.eps',trial_type));
delete(gcf);
%

%{
% repeated-measures rank tests separately for outbound and inbound
clear;
load('Wtrack_inbound_summary');
load('Wtrack_outbound_summary');
lesion_idx = find(strcmp({inbound_summary(:).group},'lesion'));
control_idx = find(strcmp({inbound_summary(:).group},'control'));
labels = {inbound_summary(:).group}';

% inbound tests are easy: we aren't missing any data
stats.in.correct = RMranktest( ...
    [inbound_summary(:).proportion_correct]',labels);
stats.in.total = RMranktest( ...
    [inbound_summary(:).number_trials]',labels);
stats.in.speed = RMranktest( ...
    [inbound_summary(:).mean_running_speed]',labels);
stats.in.duration = RMranktest( ...
    [inbound_summary(:).median_visit_duration]',labels);

% outbound, we can test total number of trials performed, but we have
% constant-size groups only for days 6-10
stats.out.total = RMranktest( ...
    [outbound_summary(:).number_trials]',labels);

correct_out = [outbound_summary(:).proportion_correct]';
correct_out = correct_out(:,6:10);
stats.out.correct = RMranktest(correct_out,labels);

speed_out = [outbound_summary(:).mean_running_speed]';
speed_out = speed_out(:,6:10);
stats.out.speed = RMranktest(speed_out,labels);

duration_out = [outbound_summary(:).median_visit_duration]';
duration_out = duration_out(:,6:10);
stats.out.duration = RMranktest(duration_out,labels);
%}

%{
% summary median+range plots for linear track
clear;
load('lineartrack_summary');
lesion_idx = find(strcmp({lineartrack_summary(:).group},'lesion'));
control_idx = find(strcmp({lineartrack_summary(:).group},'control'));

FONTSIZE = 14;
axes_style.Color = 'none';
axes_style.XLim = [0.5 2.5];
axes_style.XTick = 1:2;
axes_style.XTickLabel = {'Pre','Post'};
axes_style.FontSize = FONTSIZE;
axes_style.FontName = 'Arial';
axes_style.TickDir = 'out';
axes_style.TickLength = [0.04 0.04];
axes_style.Box = 'off';
axes_style.Color = 'none';
axes_style.LineWidth = 2;

lesion_style.Color = 'r';
lesion_style.Marker = 'o';
lesion_style.MarkerSize = 6;
lesion_style.MarkerFaceColor = 'w';
lesion_style.LineStyle = 'none';
lesion_style.LineWidth = 2;
control_style.Color = 'k';
control_style.Marker = 's';
control_style.MarkerSize = 6;
control_style.MarkerFaceColor = 'k';
control_style.LineStyle = 'none';
control_style.LineWidth = 2;

xdata = (axes_style.XTick)';
lesion_xoffset = +0.1;
control_xoffset = -0.1;

correct = vertcat(lineartrack_summary(:).correct)' ./ ...
    vertcat(lineartrack_summary(:).total)';
total = vertcat(lineartrack_summary(:).total)';
speed = vertcat(lineartrack_summary(:).overallspeed)';
duration = vertcat(lineartrack_summary(:).visitduration)';

ax_1 = axes(axes_style,'YLim',[0 1],'YTick',0:0.2:1, ...
  'Position',[0.1 0.6 0.2 0.4]);
set(get(ax_1,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Proportion correct trials');
line(xdata,median(correct(:,lesion_idx),2),'LineStyle','-','LineWidth',6, ...
    'Color',lesion_style.Color);
line(xdata,median(correct(:,control_idx),2),'LineStyle','-','LineWidth',6, ...
    'Color',control_style.Color);
line(xdata + lesion_xoffset,correct(:,lesion_idx),lesion_style);
line(xdata + control_xoffset,correct(:,control_idx),control_style);

ax_2 = axes(axes_style,'YLim',[0 250],'YTick',0:50:250, ...
    'Position',[0.5 0.6 0.2 0.4]);
set(get(ax_2,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Number of trials performed');
line(xdata,median(total(:,lesion_idx),2),'LineStyle','-','LineWidth',6, ...
    'Color',lesion_style.Color);
line(xdata,median(total(:,control_idx),2),'LineStyle','-','LineWidth',6, ...
    'Color',control_style.Color);
line(xdata + lesion_xoffset,total(:,lesion_idx),lesion_style);
line(xdata + control_xoffset,total(:,control_idx),control_style);

ax_3 = axes(axes_style,'YLim',[0 80],'YTick',0:20:80, ...
    'Position',[0.1 0.1 0.2 0.4]);
set(get(ax_3,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Running speed (cm/s)');
line(xdata,median(speed(:,lesion_idx),2),'LineStyle','-','LineWidth',6, ...
    'Color',lesion_style.Color);
line(xdata,median(speed(:,control_idx),2),'LineStyle','-','LineWidth',6, ...
    'Color',control_style.Color);
line(xdata + lesion_xoffset,speed(:,lesion_idx),lesion_style);
line(xdata + control_xoffset,speed(:,control_idx),control_style);

ax_4 = axes(axes_style,'YLim',[0 15],'YTick',0:3:15, ...
    'Position',[0.5 0.1 0.2 0.4]);
set(get(ax_4,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Average duration of food-well visit (s)');
line(xdata,median(duration(:,lesion_idx),2),'LineStyle','-','LineWidth',6, ...
    'Color',lesion_style.Color);
line(xdata,median(duration(:,control_idx),2),'LineStyle','-','LineWidth',6, ...
    'Color',control_style.Color);
line(xdata + lesion_xoffset,duration(:,lesion_idx),lesion_style);
line(xdata + control_xoffset,duration(:,control_idx),control_style);

pause;
print(gcf,'-depsc','-r600','/home/smkim/lesion_manuscript/revisions/linear_track.eps');
delete(gcf);
%}

%{
clear;
load('/data15/smkim/lesion_study/histology/estimated_volumes.mat');

control_idx = find(strcmp({volume(:).group},'control'));
lesion_idx = find(strcmp({volume(:).group},'lesion'));

FONTSIZE = 14;
axes_style.Color = 'none';
axes_style.FontSize = FONTSIZE;
axes_style.FontName = 'Arial';
axes_style.XLim = [0 1.1];
axes_style.XTick = 0:0.2:1;
axes_style.YLim = [0 1.1];
axes_style.YTick = 0:0.2:1;
axes_style.DataAspectRatio = [1 1 1];
axes_style.TickDir = 'out';
axes_style.TickLength = [0.04 0.04];
axes_style.Box = 'off';
axes_style.Color = 'none';
axes_style.LineWidth = 2;

lesion_style.Color = 'r';
lesion_style.Marker = 'o';
lesion_style.MarkerSize = 6;
lesion_style.MarkerFaceColor = 'w';
lesion_style.LineStyle = 'none';
lesion_style.LineWidth = 2;
control_style.Color = 'k';
control_style.Marker = 's';
control_style.MarkerSize = 6;
control_style.MarkerFaceColor = 'k';
control_style.LineStyle = 'none';
control_style.LineWidth = 2;


hipp_volume = ([volume(:).CA_DG_left] + [volume(:).CA_DG_right])/2;
retrohipp_volume = ([volume(:).EC_SUB_left] + [volume(:).EC_SUB_right])/2;

hipp_volume_control_mean = mean(hipp_volume(control_idx));
retrohipp_volume_control_mean = mean(retrohipp_volume(control_idx));

ax = axes('Position',[0.2 0.2 0.4 0.4],axes_style);
set(get(ax,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Normalized volume of dentate gyrus, CA fields and fimbria');
set(get(ax,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Normalized volume of retrohippocampal cortex');
line(hipp_volume(control_idx)/hipp_volume_control_mean, ...
    retrohipp_volume(control_idx)/retrohipp_volume_control_mean, ...
    control_style);
line(hipp_volume(lesion_idx)/hipp_volume_control_mean, ...
    retrohipp_volume(lesion_idx)/retrohipp_volume_control_mean, ...
    lesion_style);
line(axes_style.XLim,[1 1],'LineStyle','--','Color','k','LineWidth',1);
line([1 1],axes_style.YLim,'LineStyle','--','Color','k','LineWidth',1);

pause;
print(gcf,'-dpsc','-r600', ...
    '/home/smkim/lesion_manuscript/revisions/lesion_volumes.ps');
delete(gcf);
%}

