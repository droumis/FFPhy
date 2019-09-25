
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot example learning curves per trial, separately for inbound and outbound
clear('axes_style');
load('behavperform_2');
lesion_idx = find(strcmp({behavperform(:).group},'lesion'));
control_idx = find(strcmp({behavperform(:).group},'control'));
% maximum number of trials on either inbound or outbound for any subject
max_trials_anywhere = 1100;

FONTSIZE = 14;

for s = 1:length(behavperform)
  x = behavperform(s);
  criterion = x.criterion;
  max_in = length(x.inreward);
  max_out = length(x.outreward);
  if max_in > max_out
    max_trials = max_in;
  else
    max_trials = max_out;
  end
  if strcmp(x.group,'control')
    subject_label = sprintf('Control subject #%d',find(s==control_idx));
  elseif strcmp(x.group,'lesion')
    subject_label = sprintf('Hippocampal lesion subject #%d',find(s==lesion_idx));
  else
    error('i made a mistake');
  end

  left_margin = 10;
  bottom_margin = 0.05;
  top_margin = 0.1;

  in_axes_style.XLim = [-left_margin 100*ceil(max_in/100)];
  in_axes_style.YLim = [0-bottom_margin 1];
  in_axes_style.FontSize = FONTSIZE;
  in_axes_style.FontName = 'Arial';
  in_axes_style.XTick = [0:100:in_axes_style.XLim(2)];
  in_axes_style.YTick = 0:0.2:1;
  in_axes_style.LineWidth = 2;
  in_axes_style.TickDir = 'out';
  in_axes_style.TickLength = [0.015 0.015]* ...
      (max_trials_anywhere+left_margin)/diff(in_axes_style.XLim);
  in_axes_style.Box = 'off';
  in_axes_style.Layer = 'top';
  in_axes_style.Color = 'none';
  in_axes_style.Position = [0.15 ...
      0.57 ...
      0.75*diff(in_axes_style.XLim)/(max_trials_anywhere+left_margin) ...
      0.28];

  in_day_axes_style.Color = 'none';
  in_day_axes_style.XLim = in_axes_style.XLim;
  in_day_axes_style.YLim = [0-bottom_margin 1+top_margin];
  in_day_axes_style.Visible = 'off';
  in_day_axes_style.Position = [in_axes_style.Position(1) ...
      in_axes_style.Position(2) ...
      in_axes_style.Position(3) ...
      (top_margin + bottom_margin + 1)/(bottom_margin + 1)*in_axes_style.Position(4) ];

  out_axes_style.XLim = [-left_margin 100*ceil(max_out/100)];
  out_axes_style.YLim = [0-bottom_margin 1];
  out_axes_style.FontSize = FONTSIZE;
  out_axes_style.FontName = in_axes_style.FontName;
  out_axes_style.XTick = [0:100:out_axes_style.XLim(2)];
  out_axes_style.YTick = in_axes_style.YTick;
  out_axes_style.LineWidth = in_axes_style.LineWidth;
  out_axes_style.TickDir = in_axes_style.TickDir;
  out_axes_style.TickLength = in_axes_style.TickLength;
  out_axes_style.Box = in_axes_style.Box;
  out_axes_style.Layer = in_axes_style.Layer;
  out_axes_style.Color = 'none';
  out_axes_style.Position = [in_axes_style.Position(1) ...
      0.10 ...
      0.75*diff(out_axes_style.XLim)/(max_trials_anywhere+left_margin) ...
      in_axes_style.Position(4)];
 
  out_day_axes_style.Color = 'none'; 
  out_day_axes_style.XLim = out_axes_style.XLim;
  out_day_axes_style.YLim = in_day_axes_style.YLim;
  out_day_axes_style.Visible = 'off';
  out_day_axes_style.Position = [out_axes_style.Position(1) ...
      out_axes_style.Position(2) ...
      out_axes_style.Position(3) ...
      (top_margin + bottom_margin + 1)/(bottom_margin + 1)*out_axes_style.Position(4) ];

  framing_axes_style.Position = [0 0 1 1];
  framing_axes_style.XLim = [0 1];
  framing_axes_style.YLim = [0 1];
  framing_axes_style.Visible = 'off';
  framing_axes_style.Clipping = 'off';
  framing_axes_style.Color = 'none';

  point_style.Color = 'k';
  point_style.LineStyle = '-';
  point_style.LineWidth = 1;
  point_style.Marker = 'o';
  point_style.MarkerFaceColor = 'k';
  point_style.MarkerSize = 2;

  CI_style.FaceColor = [0.6 0.6 0.6];
  CI_style.LineStyle = 'none';

  lightgreen = [210 255 210]/255;
  darkgreen = [0 0.8 0];
  lightblue = [200 200 255]/255;
  darkblue = [0 0 0.8];
    
  criterion_color = [1 0 0];

  framing_ax = axes(framing_axes_style);
  axes(framing_ax);
  text(0,1,subject_label,'Color','k','FontSize',FONTSIZE, ...
      'FontWeight','bold', ...
      'HorizontalAlignment','Left','VerticalAlignment','top');
  text(0.05, 0.5*(in_day_axes_style.Position(2) + ...
      (out_day_axes_style.Position(2)+out_day_axes_style.Position(4))), ...
      'Estimated probability correct', ...
      'FontSize',FONTSIZE,'Color','k','Rotation',90, ...
      'HorizontalAlignment','center','VerticalAlignment','bottom');

  % draw line and confidence interval for inbound, with days indicated in
  % background
  in_ax = axes(in_axes_style);
  in_ax_right_side = axes(in_axes_style,'YAxisLocation','right');
  in_day_ax = axes(in_day_axes_style);
  axes(in_day_ax);
  text(-left_margin, 1+top_margin,'Day', ...
      'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
      'HorizontalAlignment','right','VerticalAlignment','bottom');
  for k = 1:10
    trials_range = x.dayintrials((-1+2*k):(2*k),1:2);
    if all(isnan(trials_range(:)))
      % do nothing
    else
      if (min(trials_range(:)) == 1)
        start_x = -0.5;
      else
        start_x = min(trials_range(:)) - 0.5;
      end
      end_x = max(trials_range(:)) + 0.5;
      if mod(k,2)
        patch_color = lightblue;
        label_color = darkblue;
      else
        patch_color = lightgreen;
        label_color = darkgreen;
      end
      if (k == x.ld_in)
        label_color = criterion_color;
      end
      axes(in_day_ax);
      patch([start_x start_x (start_x+end_x)/2 end_x end_x], ...
          [0-bottom_margin 1 1+top_margin 1 0-bottom_margin], ...
          patch_color);
      text((start_x + end_x)/2,1+top_margin,num2str(k), ...
          'Color',label_color,'FontSize',FONTSIZE,'FontName','Arial', ...
          'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
  end
  set(get(in_ax,'YLabel'),'FontSize',FONTSIZE,'Color','k', ...
      'String','Inbound');
  set(get(in_ax,'XLabel'),'String', ...
      'Cumulative count of inbound trials', ...
      'FontSize',FONTSIZE,'FontName','Arial','Color','k');
  axes(in_ax);
  patch('Vertices',[[(0:max_in)'; (max_in:-1:0)'], ...
      [x.inprobcorrect(:,2); flipud(x.inprobcorrect(:,3))]], ...
      'Faces',1:(2*max_in),CI_style);
  line(0:max_in,x.inprobcorrect(:,1),point_style);

  % mark learning trial
  if ~isnan(x.ld_in)
    axes(in_ax);
    line([x.lt_in x.lt_in],x.inprobcorrect(x.lt_in,2:3), ...
        'Color',criterion_color,'LineStyle','-','LineWidth',1,'Marker','none');
    line(x.lt_in,x.inprobcorrect(x.lt_in,1),point_style, ...
        'Color',criterion_color,'MarkerFaceColor',criterion_color);
    line([x.lt_in x.lt_in],[-bottom_margin 0], ...
        'Color',criterion_color,'LineStyle','-','LineWidth',1);
    text(x.lt_in,0,num2str(x.lt_in), ...
        'Color',criterion_color,'FontSize',FONTSIZE,'FontName','Arial', ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
  end

  % draw line and confidence interval for outbound, with days indicated in
  % background
  out_ax = axes(out_axes_style);
  out_ax_right_side = axes(out_axes_style,'YAxisLocation','right');
  out_day_ax = axes(out_day_axes_style);
  axes(out_day_ax);
  text(-left_margin, 1+top_margin,'Day', ...
      'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
      'HorizontalAlignment','right','VerticalAlignment','bottom');
  for k = 1:10
    trials_range = x.dayouttrials((-1+2*k):(2*k),1:2);
    if all(isnan(trials_range(:)))
      % do nothing
    else
      if (min(trials_range(:)) == 1)
        start_x = -0.5;
      else
        start_x = min(trials_range(:)) - 0.5;
      end
      end_x = max(trials_range(:)) + 0.5;
      if mod(k,2)
        patch_color = lightblue;
        label_color = darkblue;
      else
        patch_color = lightgreen;
        label_color = darkgreen;
      end
      if (k == x.ld_out)
        label_color = criterion_color;
      end
      axes(out_day_ax);
      patch([start_x start_x (start_x+end_x)/2 end_x end_x], ...
          [0-bottom_margin 1 1+top_margin 1 0-bottom_margin], ...
          patch_color);
      
      text((start_x + end_x)/2,1+top_margin,num2str(k), ...
          'Color',label_color,'FontSize',FONTSIZE,'FontName','Arial', ...
          'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
  end
  set(get(out_ax,'YLabel'),'FontSize',FONTSIZE,'Color','k', ...
      'String','Outbound');
  set(get(out_ax,'XLabel'),'String', ...
      'Cumulative count of outbound trials', ...
      'FontSize',FONTSIZE,'FontName','Arial','Color','k');
  axes(out_ax);
  patch('Vertices',[[(0:max_out)'; (max_out:-1:0)'], ...
      [x.outprobcorrect(:,2); flipud(x.outprobcorrect(:,3))]], ...
      'Faces',1:(2*max_out),CI_style);
  line(0:max_out,x.outprobcorrect(:,1),point_style);

  if ~isnan(x.ld_out)
    axes(out_ax);
    line([x.lt_out x.lt_out],x.outprobcorrect(x.lt_out,2:3), ...
        'Color',criterion_color,'LineStyle','-','LineWidth',1,'Marker','none');
    line(x.lt_out,x.outprobcorrect(x.lt_out,1),point_style, ...
        'Color',criterion_color,'MarkerFaceColor',criterion_color);
    line([x.lt_out x.lt_out],[-bottom_margin 0], ...
        'Color',criterion_color,'LineStyle','-','LineWidth',1);
    text(x.lt_out,0,num2str(x.lt_out), ...
        'Color',criterion_color,'FontSize',FONTSIZE,'FontName','Arial', ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
  end

  % dashed lines to indicate background probability
  axes(in_ax);
  line(get(gca,'XLim'),[criterion criterion],'Color','k','LineWidth',1,'LineStyle','--');
  axes(out_ax);
  line(get(gca,'XLim'),[criterion criterion],'Color','k','LineWidth',1,'LineStyle','--');

  if strcmp(x.group,'control')
    filename = sprintf( ...
        '/home/smkim/lesion_manuscript/revisions/control%d_criterion(%f).eps', ...
        find(s==control_idx),x.criterion);
  elseif strcmp(x.group,'lesion')
    filename = sprintf( ...
        '/home/smkim/lesion_manuscript/revisions/lesion%d_criterion(%f).eps', ...
        find(s==lesion_idx),x.criterion);
  else
    error('i made a mistake');
  end
  print(gcf,'-depsc','-r600',filename);
  pause;
  delete(gcf);
end
%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot example MA curves, separately for inbound and outbound
clear('axes_style');
load('behavperform_2');
lesion_idx = find(strcmp({behavperform(:).group},'lesion'));
control_idx = find(strcmp({behavperform(:).group},'control'));
% maximum number of trials on either inbound or outbound for any subject
max_trials_anywhere = 1100;

MA_order = 10;

FONTSIZE = 14;

for s = 1:length(behavperform)
  x = behavperform(s);
  criterion = x.criterion;
  max_in = length(x.inreward);
  max_out = length(x.outreward);
  if max_in > max_out
    max_trials = max_in;
  else
    max_trials = max_out;
  end
  if strcmp(x.group,'control')
    subject_label = sprintf('Control subject #%d',find(s==control_idx));
  elseif strcmp(x.group,'lesion')
    subject_label = sprintf('Hippocampal lesion subject #%d',find(s==lesion_idx));
  else
    error('i made a mistake');
  end

  left_margin = 10;
  bottom_margin = 0.05;
  top_margin = 0.1;

  in_axes_style.XLim = [-left_margin 100*ceil(max_in/100)];
  in_axes_style.YLim = [0-bottom_margin 1];
  in_axes_style.FontSize = FONTSIZE;
  in_axes_style.FontName = 'Arial';
  in_axes_style.XTick = [0:100:in_axes_style.XLim(2)];
  in_axes_style.YTick = 0:0.2:1;
  in_axes_style.LineWidth = 2;
  in_axes_style.TickDir = 'out';
  in_axes_style.TickLength = [0.015 0.015]* ...
      (max_trials_anywhere+left_margin)/diff(in_axes_style.XLim);
  in_axes_style.Box = 'off';
  in_axes_style.Layer = 'top';
  in_axes_style.Color = 'none';
  in_axes_style.Position = [0.15 ...
      0.57 ...
      0.75*diff(in_axes_style.XLim)/(max_trials_anywhere+left_margin) ...
      0.28];

  in_day_axes_style.Color = 'none';
  in_day_axes_style.XLim = in_axes_style.XLim;
  in_day_axes_style.YLim = [0-bottom_margin 1+top_margin];
  in_day_axes_style.Visible = 'off';
  in_day_axes_style.Position = [in_axes_style.Position(1) ...
      in_axes_style.Position(2) ...
      in_axes_style.Position(3) ...
      (top_margin + bottom_margin + 1)/(bottom_margin + 1)*in_axes_style.Position(4) ];

  out_axes_style.XLim = [-left_margin 100*ceil(max_out/100)];
  out_axes_style.YLim = [0-bottom_margin 1];
  out_axes_style.FontSize = FONTSIZE;
  out_axes_style.FontName = in_axes_style.FontName;
  out_axes_style.XTick = [0:100:out_axes_style.XLim(2)];
  out_axes_style.YTick = in_axes_style.YTick;
  out_axes_style.LineWidth = in_axes_style.LineWidth;
  out_axes_style.TickDir = in_axes_style.TickDir;
  out_axes_style.TickLength = in_axes_style.TickLength;
  out_axes_style.Box = in_axes_style.Box;
  out_axes_style.Layer = in_axes_style.Layer;
  out_axes_style.Color = 'none';
  out_axes_style.Position = [in_axes_style.Position(1) ...
      0.10 ...
      0.75*diff(out_axes_style.XLim)/(max_trials_anywhere+left_margin) ...
      in_axes_style.Position(4)];
 
  out_day_axes_style.Color = 'none'; 
  out_day_axes_style.XLim = out_axes_style.XLim;
  out_day_axes_style.YLim = in_day_axes_style.YLim;
  out_day_axes_style.Visible = 'off';
  out_day_axes_style.Position = [out_axes_style.Position(1) ...
      out_axes_style.Position(2) ...
      out_axes_style.Position(3) ...
      (top_margin + bottom_margin + 1)/(bottom_margin + 1)*out_axes_style.Position(4) ];

  framing_axes_style.Position = [0 0 1 1];
  framing_axes_style.XLim = [0 1];
  framing_axes_style.YLim = [0 1];
  framing_axes_style.Visible = 'off';
  framing_axes_style.Clipping = 'off';
  framing_axes_style.Color = 'none';

  point_style.Color = 'k';
  point_style.LineStyle = '-';
  point_style.LineWidth = 1;
  point_style.Marker = 'o';
  point_style.MarkerFaceColor = 'k';
  point_style.MarkerSize = 2;

  CI_style.FaceColor = [0.6 0.6 0.6];
  CI_style.LineStyle = 'none';

  lightgreen = [210 255 210]/255;
  darkgreen = [0 0.8 0];
  lightblue = [200 200 255]/255;
  darkblue = [0 0 0.8];
    
  criterion_color = [1 0 0];

  framing_ax = axes(framing_axes_style);
  axes(framing_ax);
  text(0,1,subject_label,'Color','k','FontSize',FONTSIZE, ...
      'FontWeight','bold', ...
      'HorizontalAlignment','Left','VerticalAlignment','top');
  text(0.05, 0.5*(in_day_axes_style.Position(2) + ...
      (out_day_axes_style.Position(2)+out_day_axes_style.Position(4))), ...
      sprintf('Proportion correct in %d-trial moving window',MA_order), ...
      'FontSize',FONTSIZE,'Color','k','Rotation',90, ...
      'HorizontalAlignment','center','VerticalAlignment','bottom');

  % draw moving average trendline for inbound, with days indicated in
  % background
  in_ax = axes(in_axes_style);
  in_ax_right_side = axes(in_axes_style,'YAxisLocation','right');
  in_day_ax = axes(in_day_axes_style);
  axes(in_day_ax);
  text(-left_margin, 1+top_margin,'Day', ...
      'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
      'HorizontalAlignment','right','VerticalAlignment','bottom');
  for k = 1:10
    trials_range = x.dayintrials((-1+2*k):(2*k),1:2);
    if all(isnan(trials_range(:)))
      % do nothing
    else
      if (min(trials_range(:)) == 1)
        start_x = -0.5;
      else
        start_x = min(trials_range(:)) - 0.5;
      end
      end_x = max(trials_range(:)) + 0.5;
      if mod(k,2)
        patch_color = lightblue;
        label_color = darkblue;
      else
        patch_color = lightgreen;
        label_color = darkgreen;
      end
      axes(in_day_ax);
      patch([start_x start_x (start_x+end_x)/2 end_x end_x], ...
          [0-bottom_margin 1 1+top_margin 1 0-bottom_margin], ...
          patch_color);
      text((start_x + end_x)/2,1+top_margin,num2str(k), ...
          'Color',label_color,'FontSize',FONTSIZE,'FontName','Arial', ...
          'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
  end
  set(get(in_ax,'YLabel'),'FontSize',FONTSIZE,'Color','k', ...
      'String','Inbound');
  set(get(in_ax,'XLabel'),'String', ...
      'Cumulative count of inbound trials', ...
      'FontSize',FONTSIZE,'FontName','Arial','Color','k');
  axes(in_day_ax);
  inbound_ma = convn(x.inreward,ones(MA_order,1)/MA_order,'valid');
  line(MA_order:max_in,inbound_ma,point_style);
  axes(in_ax);

  % draw moving average trendline for outbound, with days indicated in
  % background
  out_ax = axes(out_axes_style);
  out_ax_right_side = axes(out_axes_style,'YAxisLocation','right');
  out_day_ax = axes(out_day_axes_style);
  axes(out_day_ax);
  text(-left_margin, 1+top_margin,'Day', ...
      'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
      'HorizontalAlignment','right','VerticalAlignment','bottom');
  for k = 1:10
    trials_range = x.dayouttrials((-1+2*k):(2*k),1:2);
    if all(isnan(trials_range(:)))
      % do nothing
    else
      if (min(trials_range(:)) == 1)
        start_x = -0.5;
      else
        start_x = min(trials_range(:)) - 0.5;
      end
      end_x = max(trials_range(:)) + 0.5;
      if mod(k,2)
        patch_color = lightblue;
        label_color = darkblue;
      else
        patch_color = lightgreen;
        label_color = darkgreen;
      end
      axes(out_day_ax);
      patch([start_x start_x (start_x+end_x)/2 end_x end_x], ...
          [0-bottom_margin 1 1+top_margin 1 0-bottom_margin], ...
          patch_color);
      
      text((start_x + end_x)/2,1+top_margin,num2str(k), ...
          'Color',label_color,'FontSize',FONTSIZE,'FontName','Arial', ...
          'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
  end
  set(get(out_ax,'YLabel'),'FontSize',FONTSIZE,'Color','k', ...
      'String','Outbound');
  set(get(out_ax,'XLabel'),'String', ...
      'Cumulative count of outbound trials', ...
      'FontSize',FONTSIZE,'FontName','Arial','Color','k');
  axes(out_day_ax);
  outbound_ma = convn(x.outreward,ones(MA_order,1)/MA_order,'valid');
  line(MA_order:max_out,outbound_ma,point_style);
  axes(out_ax);

  if strcmp(x.group,'control')
    filename = sprintf( ...
        '/home/smkim/lesion_manuscript/revisions/control%d_MA%d.eps', ...
        find(s==control_idx),MA_order);
  elseif strcmp(x.group,'lesion')
    filename = sprintf( ...
        '/home/smkim/lesion_manuscript/revisions/lesion%d_MA%d.eps', ...
        find(s==lesion_idx),MA_order);
  else
    error('i made a mistake');
  end
  print(gcf,'-depsc','-r600',filename);
  pause;
  delete(gcf);
end
%}

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tabulate various things about the learning curves, 
% separately for inbound and outbound, for every subject, 
% with between-group p-values:
%   total number of trials
%   trials to reach criterion
%   day on which learning criterion was achieved
%   average estimated probability correct on day 10
%   
%   
clear;
load('behavperform_2');
lesion_idx = find(strcmp({behavperform(:).group},'lesion'));
control_idx = find(strcmp({behavperform(:).group},'control'));
lastday = 10;
halfwidth = 5;

for s = 1:length(behavperform)
  x = behavperform(s);

  if strcmp(x.group,'control')
    learning_stats(s).subject_label = ...
        sprintf('Control subject #%d',find(s==control_idx));
  elseif strcmp(x.group,'lesion')
    learning_stats(s).subject_label =  ...
        sprintf('Hippocampal lesion subject #%d',find(s==lesion_idx));
  else
    error('i made a mistake');
  end
  learning_stats(s).group = x.group;
  learning_stats(s).criterion = x.criterion;
  learning_stats(s).trialcount_in = length(x.inreward);
  learning_stats(s).ld_in = x.ld_in;
  learning_stats(s).lastday_in = mean(x.inprobcorrect( ...
    (x.dayintrials(2*lastday-1,1)):(x.dayintrials(2*lastday,2))',1));
  if ~isnan(x.ld_in)
    learning_stats(s).lt_in = x.lt_in;
  else
    lastbelowcriterion = find(x.inprobcorrect(:,2) <= x.criterion,1,'last');
    if (lastbelowcriterion <= learning_stats(s).trialcount_in)
      learning_stats(s).lt_in = lastbelowcriterion;
    else
      learning_stats(s).lt_in = learning_stats(s).trialcount_in + 1;
    end
  end

  learning_stats(s).trialcount_out = length(x.outreward);
  learning_stats(s).ld_out = x.ld_out;
  learning_stats(s).lastday_out = mean(x.outprobcorrect( ...
    (x.dayouttrials(2*lastday-1,1)):(x.dayouttrials(2*lastday,2))',1));
  if ~isnan(x.ld_out)
    learning_stats(s).lt_out = x.lt_out;
  else
    lastbelowcriterion = find(x.outprobcorrect(:,2) <= x.criterion,1,'last');
    if (lastbelowcriterion <= learning_stats(s).trialcount_in)
      learning_stats(s).lt_out = lastbelowcriterion;
    else
      learning_stats(s).lt_out = learning_stats(s).trialcount_in + 1;
    end
  end
end
save('Wtrack_inbound_outbound_stats','learning_stats');

%
FONTSIZE = 14;

lesion_style.Color = 'r';
lesion_style.LineStyle = 'none';
lesion_style.Marker = 'o';
lesion_style.LineWidth = 2;
lesion_style.MarkerSize = 6;
lesion_style.MarkerFaceColor = 'none';
control_style.Color = 'k';
control_style.LineStyle = 'none';
control_style.Marker = 's';
control_style.LineWidth = 2;
control_style.MarkerSize = 6;
control_style.MarkerFaceColor = 'k';

% plot days to reach learning criterion on inbound and outbound components
axes_style.YLim = [0.5 2.5];
axes_style.XLim = [0.5 11];
axes_style.XTick = [1:8, 8.5, 9.5, 10.5];
axes_style.XTickLabel = horzcat(arrayfun(@num2str,1:8, ...
    'UniformOutput',false),{'','Did not reach criterion',''});
axes_style.YTick = 1:2;
axes_style.YTickLabel = {'Hippocampal lesion','Control'};
axes_style.FontSize = FONTSIZE;
axes_style.FontName = 'Arial';
axes_style.LineWidth = 2;
axes_style.TickDir = 'out';
axes_style.TickLength = [0.02 0.02];
axes_style.Box = 'off';
axes_style.Layer = 'top';
axes_style.Color = 'none';
axes_style.YDir = 'normal';

ld_out = [learning_stats(:).ld_out];
ld_out(isnan(ld_out)) = 10.5;
ld_in = [learning_stats(:).ld_in];
ld_in(isnan(ld_in)) = 10.5;
ld_in_ydata(control_idx) = 2*ones(numel(control_idx),1);
ld_in_ydata(lesion_idx) = 1*ones(numel(lesion_idx),1);
ld_out_ydata(control_idx) = 2*ones(numel(control_idx),1);
ld_out_ydata(lesion_idx) = 1*ones(numel(lesion_idx),1);
jitter = 0.25;
dayrange = axes_style.XTick;
for i = 1:numel(dayrange)
    idx = control_idx(find(ld_in(control_idx) == dayrange(i)));
    ld_in_ydata(idx) = ld_in_ydata(idx) + jitter*((1:numel(idx))-(1+numel(idx))/2);
    idx = lesion_idx(find(ld_in(lesion_idx) == dayrange(i)));
    ld_in_ydata(idx) = ld_in_ydata(idx) + jitter*((1:numel(idx))-(1+numel(idx))/2);
    idx = control_idx(find(ld_out(control_idx) == dayrange(i)));
    ld_out_ydata(idx) = ld_out_ydata(idx) + jitter*((1:numel(idx))-(1+numel(idx))/2);
    idx = lesion_idx(find(ld_out(lesion_idx) == dayrange(i)));
    ld_out_ydata(idx) = ld_out_ydata(idx) + jitter*((1:numel(idx))-(1+numel(idx))/2);
end
figure(1);
a1 = axes('Position',[0.3 0.5 0.5 0.2],axes_style);
set(get(a1,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Days to reach learning criterion');
text(6.5,2.5,'Inbound','FontSize',FONTSIZE,'FontName','Arial', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center');
legend_line_handles(1,1) = line(ld_in(control_idx), ...
    ld_in_ydata(control_idx),control_style);
legend_line_handles(2,1) = line(ld_in(lesion_idx), ...
    ld_in_ydata(lesion_idx),lesion_style);
% legend
[legend_h,object_h,plot_h,text_strings] = ...
    legend(legend_line_handles(:,1), ...
    sprintf('Control (%d rats)',nnz(control_idx)), ...
    sprintf('Hippocampal lesion (%d rats)',nnz(lesion_idx)), ...
    'Location','none');
set(legend_h,'Position',[0.1 0.78 0.2 0.1]);
set(object_h(1),'FontSize',FONTSIZE,'FontName','Arial');
set(object_h(2),'FontSize',FONTSIZE,'FontName','Arial');
legend(a1,'boxoff');

a2 = axes('Position',[0.3 0.1 0.5 0.2],axes_style);
set(get(a2,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Days to reach learning criterion');
text(6.5,2.5,'Outbound','FontSize',FONTSIZE,'FontName','Arial', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center');
line(ld_out(control_idx),ld_out_ydata(control_idx),control_style);
line(ld_out(lesion_idx),ld_out_ydata(lesion_idx),lesion_style);
% ranksum test on days to reach learning criterion
disp(sprintf('days to criterion, ranksum test:  inbound p %f  outbound p %f', ...
    [ ranksum(ld_in(lesion_idx),ld_in(control_idx)) ...
    ranksum(ld_out(lesion_idx),ld_out(control_idx)) ]));
print(gcf,'-depsc','-r600','/home/smkim/lesion_manuscript/revisions/days_to_criterion.eps');

% plot "asymptotic" performance
axes_style.YLim = [0.5 2.5];
axes_style.XLim = [0 1];
axes_style.XTick = 0:0.1:1;
axes_style.XTickLabel = arrayfun(@num2str,axes_style.XTick, ...
    'UniformOutput',false);
axes_style.YTick = 1:2;
axes_style.YTickLabel = {'Hippocampal lesion','Control'};
axes_style.FontSize = FONTSIZE;
axes_style.FontName = 'Arial';
axes_style.LineWidth = 2;
axes_style.TickDir = 'out';
axes_style.TickLength = [0.02 0.02];
axes_style.Box = 'off';
axes_style.Layer = 'top';
axes_style.Color = 'none';
axes_style.YDir = 'normal';

lastday_out = [learning_stats(:).lastday_out];
lastday_in = [learning_stats(:).lastday_in];
lastday_in_ydata(control_idx) = 2*ones(numel(control_idx),1);
lastday_in_ydata(lesion_idx) = 1*ones(numel(lesion_idx),1);
lastday_out_ydata(control_idx) = 2*ones(numel(control_idx),1);
lastday_out_ydata(lesion_idx) = 1*ones(numel(lesion_idx),1);
jitter = 0.10;
lastday_out_ydata = lastday_out_ydata + jitter*randn(size(lastday_out_ydata));
lastday_in_ydata = lastday_in_ydata + jitter*randn(size(lastday_in_ydata));
figure(2);
a1 = axes('Position',[0.3 0.5 0.5 0.2],axes_style);
set(get(a1,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String', ...
    sprintf('Mean estimated probability correct on day %d',lastday));
text(0.5,2.5,'Inbound','FontSize',FONTSIZE,'FontName','Arial', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center');
legend_line_handles(1,2) = line(lastday_in(control_idx), ...
    lastday_in_ydata(control_idx),control_style);
legend_line_handles(2,2) = line(lastday_in(lesion_idx), ...
    lastday_in_ydata(lesion_idx),lesion_style);
% legend
[legend_h,object_h,plot_h,text_strings] = ...
    legend(legend_line_handles(:,2), ...
    sprintf('Control (%d rats)',nnz(control_idx)), ...
    sprintf('Hippocampal lesion (%d rats)',nnz(lesion_idx)), ...
    'Location','none');
set(legend_h,'Position',[0.1 0.78 0.2 0.1]);
set(object_h(1),'FontSize',FONTSIZE,'FontName','Arial');
set(object_h(2),'FontSize',FONTSIZE,'FontName','Arial');
legend(a1,'boxoff');

a2 = axes('Position',[0.3 0.1 0.5 0.2],axes_style);
set(get(a2,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String', ...
    sprintf('Mean estimated probability correct on day %d',lastday));
text(0.5,2.5,'Outbound','FontSize',FONTSIZE,'FontName','Arial', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center');
line(lastday_out(control_idx),lastday_out_ydata(control_idx),control_style);
line(lastday_out(lesion_idx),lastday_out_ydata(lesion_idx),lesion_style);
% ranksum test on asymptotic performance
disp(sprintf('day %d performance, ranksum test:  inbound p %f  outbound p %f', ...
    [ lastday ranksum(lastday_in(lesion_idx),lastday_in(control_idx)) ...
    ranksum(lastday_out(lesion_idx),lastday_out(control_idx)) ]));
print(gcf,'-depsc','-r600','/home/smkim/lesion_manuscript/revisions/lastday_performance.eps');
%
pause;
delete(1); delete(2);
%

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
      'median_running_speed',nan(10,1), ...
      'median_visit_duration',nan(10,1));
  outbound_summary(s) = struct(...
      'subject',subject_label, ...
      'group',behavperform(s).group, ...
      'number_trials',nan(10,1), ...
      'proportion_correct',nan(10,1), ...
      'proportion_revisit_error',nan(10,1), ...
      'proportion_side_error',nan(10,1), ...
      'median_running_speed',nan(10,1), ...
      'median_visit_duration',nan(10,1));
  for d = 1:10
    inbound_trials = struct( ...
        'origin',{}, ...
        'destination',{}, ...
        'correct',{}, ...
        'running_speed',{}, ...
        'visit_duration',{});
    outbound_trials = struct( ...
        'origin',{}, ...
        'destination',{}, ...
        'correct',{}, ...
        'running_speed',{}, ...
        'visit_duration',{});
    for r = 1:2
      for j = 1:length(journeys{d}{r}.correct)
        if ~isnan(journeys{d}{r}.correct(j))
          if (journeys{d}{r}.startzone == 2)
            outbound_trials(end+1) = struct( ...
              'origin',journeys{d}{r}.startzone(j), ...
              'destination',journeys{d}{r}.endzone(j), ...
              'correct',journeys{d}{r}.correct(j), ...
              'running_speed',journeys{d}{r}.speed(j), ...
              'visit_duration',journeys{d}{r}.visitduration(j) );
          else
            inbound_trials(end+1) = struct( ...
              'origin',journeys{d}{r}.startzone(j), ...
              'destination',journeys{d}{r}.endzone(j), ...
              'correct',journeys{d}{r}.correct(j), ...
              'running_speed',journeys{d}{r}.speed(j), ...
              'visit_duration',journeys{d}{r}.visitduration(j) );
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
    inbound_summary(s).median_running_speed = median( ...
        [inbound_trials(:).running_speed]);
    inbound_summary(s).median_visit_duration = median( ...
        [inbound_trials(:).visit_duration]);

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
    outbound_summary(s).median_running_speed = median( ...
        [outbound_trials(:).running_speed]);
    outbound_summary(s).median_visit_duration = median( ...
        [outbound_trials(:).visit_duration]);
  end
end

save('Wtrack_inbound_summary.mat','inbound_summary');
save('Wtrack_outbound_summary.mat','outbound_summary');
%}

%{
% lower-triangle plot of inbound errors, with return errors versus side-to-side
% errors on the two axes
clear;
load('Wtrack_inbound_summary');
load('Wtrack_outbound_summary');
lesion_idx = find(strcmp({inbound_summary(:).group},'lesion'));
control_idx = find(strcmp({inbound_summary(:).group},'control'));

FONTSIZE = 14;
axes_style.XLim = [0 1];
axes_style.XTick = 0:0.25:1;
axes_style.YLim = [0 1];
axes_style.YTick = 0:0.25:1;
axes_style.FontSize = FONTSIZE;
axes_style.FontName = 'Arial';
axes_style.LineWidth = 2;
axes_style.TickDir = 'out';
axes_style.TickLength = [0.04 0.04];
axes_style.Box = 'off';
axes_style.Layer = 'bottom';
axes_style.Color = 'none';
axes_style.DataAspectRatio = [1 1 1];
axes_style.Visible = 'on';
axes_style.YDir = 'normal';

lesion_style.Color = 'r';
lesion_style.LineStyle = 'none';
lesion_style.Marker = 'o';
lesion_style.LineWidth = 2;
lesion_style.MarkerSize = 6;
lesion_style.MarkerFaceColor = 'none';
control_style.Color = 'k';
control_style.LineStyle = 'none';
control_style.Marker = 's';
control_style.LineWidth = 2;
control_style.MarkerSize = 6;
control_style.MarkerFaceColor = 'k';

a1 = axes(axes_style,'Position',[0.1 0.4 0.4 0.4]);
set(get(a1,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Proportion of side-to-side errors on inbound trials');
set(get(a1,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Proportion of turn-around errors on inbound trials');
day = 1;
text(0.5,1,sprintf('Day %d',day),'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center','FontSize',FONTSIZE);
% plot diagonal line
line([0 1],[1 0], ...
  'LineWidth',axes_style.LineWidth,'Color','k','LineStyle','--');
% plot points
inbound_side_error = horzcat(inbound_summary(:).proportion_side_error);
inbound_revisit_error = horzcat(inbound_summary(:).proportion_revisit_error);
legend_line_handles(1) = line(inbound_side_error(day,control_idx), ...
    inbound_revisit_error(day,control_idx), ...
    control_style);
legend_line_handles(2) = line(inbound_side_error(day,lesion_idx), ...
    inbound_revisit_error(day,lesion_idx), ...
    lesion_style);
% legend
[legend_h,object_h,plot_h,text_strings] = ...
    legend(legend_line_handles, ...
    sprintf('Control (%d rats)',nnz(control_idx)), ...
    sprintf('Hippocampal lesion (%d rats)',nnz(lesion_idx)), ...
    'Location','none');
set(legend_h,'Position',[0.1 0.83 0.2 0.1]);
set(object_h(1),'FontSize',FONTSIZE,'FontName','Arial');
set(object_h(2),'FontSize',FONTSIZE,'FontName','Arial');
legend(a1,'boxoff');
%
a2 = axes(axes_style,'Position',[0.55 0.4 0.4 0.4]);
set(get(a2,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Proportion of side-to-side errors on inbound trials');
set(get(a2,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Proportion of turn-around errors on inbound trials');
day = 2;
text(0.5,1,sprintf('Day %d',day),'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center','FontSize',FONTSIZE);
% plot diagonal line
line([0 1],[1 0], ...
  'LineWidth',axes_style.LineWidth,'Color','k','LineStyle','--');
inbound_side_error = horzcat(inbound_summary(:).proportion_side_error);
inbound_revisit_error = horzcat(inbound_summary(:).proportion_revisit_error);
line(inbound_side_error(day,control_idx),inbound_revisit_error(day,control_idx), ...
    control_style);
line(inbound_side_error(day,lesion_idx),inbound_revisit_error(day,lesion_idx), ...
    lesion_style);

pause;
print(gcf,'-depsc','-r600','/home/smkim/lesion_manuscript/revisions/inbound_errors.eps');
delete(gcf);
%}

%{
% summary median+range for W track, separated by inbound and outbound
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
axes_style.LineWidth = 2;
axes_style.TickDir = 'out';
axes_style.TickLength = [0.02 0.02];
axes_style.Box = 'off';
axes_style.Color = 'none';

lesion_style.Color = 'r';
lesion_style.LineStyle = '-';
lesion_style.Marker = 'o';
lesion_style.LineWidth = 2;
lesion_style.MarkerSize = 6;
lesion_style.MarkerFaceColor = 'w';
control_style.Color = 'k';
control_style.LineStyle = '-';
control_style.Marker = 's';
control_style.LineWidth = 2;
control_style.MarkerSize = 6;
control_style.MarkerFaceColor = 'k';

xdata = axes_style.XTick;
lesion_xoffset = +0.10;
control_xoffset = -0.10;

correct_in = [inbound_summary(:).proportion_correct];
correct_out = [outbound_summary(:).proportion_correct];

total_in = [inbound_summary(:).number_trials];
total_out = [outbound_summary(:).number_trials];

speed_in = [inbound_summary(:).median_running_speed];
speed_out = [outbound_summary(:).median_running_speed];

duration_in = [inbound_summary(:).median_visit_duration];
duration_out = [outbound_summary(:).median_visit_duration];

ax_1 = axes(axes_style,'YLim',[-0.05 1],'YTick',0:0.2:1, ...
  'Position',[0.1 0.5 0.8 0.32]);
text(mean(axes_style.XLim),1, ...
    'Inbound','VerticalAlignment','bottom', ...
    'HorizontalAlignment','center','FontSize',FONTSIZE);
set(get(ax_1,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Proportion correct trials');
set(get(ax_1,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Day');
line(repmat(xdata,[2 1]) + control_xoffset, ...
    [min(correct_out(control_idx,:)); max(correct_out(control_idx,:))], ...
    'LineStyle',control_style.LineStyle, ...
    'LineWidth',control_style.LineWidth, ...
    'Color',control_style.Color);
legend_line_handles(1) = line(xdata + control_xoffset, ...
    median(correct_out(control_idx,:),1),control_style);
line(repmat(xdata,[2 1]) + lesion_xoffset, ...
    [min(correct_out(lesion_idx,:)); max(correct_out(lesion_idx,:))], ...
    'LineStyle',lesion_style.LineStyle, ...
    'LineWidth',lesion_style.LineWidth, ...
    'Color',lesion_style.Color);
legend_line_handles(2) = line(xdata + lesion_xoffset, ...
    median(correct_out(lesion_idx,:),1),lesion_style);
%
[legend_h,object_h,plot_h,text_strings] = ...
    legend(legend_line_handles, ...
    sprintf('Control (%d rats)',nnz(control_idx)), ...
    sprintf('Hippocampal lesion (%d rats)',nnz(lesion_idx)), ...
    'Location','none');
set(legend_h,'Position',[0.1 0.9 0.2 0.1]);
set(object_h(1),'FontSize',FONTSIZE,'FontName','Arial');
set(object_h(2),'FontSize',FONTSIZE,'FontName','Arial');
legend(ax_1,'boxoff');
%

pause;
print(gcf,'-depsc','-r600','/home/smkim/lesion_manuscript/revisions/test.eps');
delete(gcf);
%}


%{
% summary median+range plots for linear track
clear;
load('lineartrack_summary');
lesion_idx = find(strcmp({lineartrack_summary(:).group},'lesion'));
control_idx = find(strcmp({lineartrack_summary(:).group},'control'));

FONTSIZE = 14;
axes_style.Color = 'none';
axes_style.XLim = [0 3];
axes_style.XTick = [0.5 2.5];
axes_style.XTickLabel = {'Pre','Post'};
axes_style.FontSize = FONTSIZE;
axes_style.FontName = 'Arial';
axes_style.LineWidth = 2;
axes_style.TickDir = 'out';
axes_style.TickLength = [0.04 0.04];
axes_style.Box = 'off';
axes_style.Color = 'none';

lesion_style.Color = 'r';
lesion_style.LineStyle = '-';
lesion_style.Marker = 'o';
lesion_style.LineWidth = 2;
lesion_style.MarkerSize = 6;
lesion_style.MarkerFaceColor = 'w';
control_style.Color = 'k';
control_style.LineStyle = '-';
control_style.Marker = 's';
control_style.LineWidth = 2;
control_style.MarkerSize = 6;
control_style.MarkerFaceColor = 'k';

xdata = axes_style.XTick;
lesion_xoffset = +0.10;
control_xoffset = -0.10;

correct = vertcat(lineartrack_summary(:).correct) ./ ...
    vertcat(lineartrack_summary(:).total);
total = vertcat(lineartrack_summary(:).total);
speed = vertcat(lineartrack_summary(:).overallspeed);
duration = vertcat(lineartrack_summary(:).visitduration);

ax_1 = axes(axes_style,'YLim',[0 1],'YTick',0:0.2:1, ...
  'Position',[0.1 0.5 0.2 0.4]);
set(get(ax_1,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Proportion correct trials');
line(repmat(xdata,[2 1]) + control_xoffset, ...
    [min(correct(control_idx,:)); max(correct(control_idx,:))], ...
    'LineStyle',control_style.LineStyle, ...
    'LineWidth',control_style.LineWidth, ...
    'Color',control_style.Color);
legend_line_handles(1) = line(xdata + control_xoffset, ...
    median(correct(control_idx,:),1),control_style);
line(repmat(xdata,[2 1]) + lesion_xoffset, ...
    [min(correct(lesion_idx,:)); max(correct(lesion_idx,:))], ...
    'LineStyle',lesion_style.LineStyle, ...
    'LineWidth',lesion_style.LineWidth, ...
    'Color',lesion_style.Color);
legend_line_handles(2) = line(xdata + lesion_xoffset, ...
    median(correct(lesion_idx,:),1),lesion_style);
%
[legend_h,object_h,plot_h,text_strings] = ...
    legend(legend_line_handles, ...
    sprintf('Control (%d rats)',nnz(control_idx)), ...
    sprintf('Hippocampal lesion (%d rats)',nnz(lesion_idx)), ...
    'Location','none');
set(legend_h,'Position',[0.1 0.9 0.2 0.1]);
set(object_h(1),'FontSize',FONTSIZE,'FontName','Arial');
set(object_h(2),'FontSize',FONTSIZE,'FontName','Arial');
legend(ax_1,'boxoff');
%

ax_2 = axes(axes_style,'YLim',[0 250],'YTick',0:50:250, ...
    'Position',[0.55 0.5 0.2 0.4]);
set(get(ax_2,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Number of trials performed');
line(repmat(xdata,[2 1]) + control_xoffset, ...
    [min(total(control_idx,:)); max(total(control_idx,:))], ...
    'LineStyle',control_style.LineStyle, ...
    'LineWidth',control_style.LineWidth, ...
    'Color',control_style.Color);
line(xdata + control_xoffset, ...
    median(total(control_idx,:),1),control_style);
line(repmat(xdata,[2 1]) + lesion_xoffset, ...
    [min(total(lesion_idx,:)); max(total(lesion_idx,:))], ...
    'LineStyle',lesion_style.LineStyle, ...
    'LineWidth',lesion_style.LineWidth, ...
    'Color',lesion_style.Color);
line(xdata + lesion_xoffset, ...
    median(total(lesion_idx,:),1),lesion_style);

ax_3 = axes(axes_style,'YLim',[0 80],'YTick',0:20:80, ...
    'Position',[0.1 0.1 0.2 0.4]);
set(get(ax_3,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Running speed (cm/s)');
line(repmat(xdata,[2 1]) + control_xoffset, ...
    [min(speed(control_idx,:)); max(speed(control_idx,:))], ...
    'LineStyle',control_style.LineStyle, ...
    'LineWidth',control_style.LineWidth, ...
    'Color',control_style.Color);
line(xdata + control_xoffset, ...
    median(speed(control_idx,:),1),control_style);
line(repmat(xdata,[2 1]) + lesion_xoffset, ...
    [min(speed(lesion_idx,:)); max(speed(lesion_idx,:))], ...
    'LineStyle',lesion_style.LineStyle, ...
    'LineWidth',lesion_style.LineWidth, ...
    'Color',lesion_style.Color);
line(xdata + lesion_xoffset, ...
    median(speed(lesion_idx,:),1),lesion_style);

ax_4 = axes(axes_style,'YLim',[0 15],'YTick',0:3:15, ...
    'Position',[0.55 0.1 0.2 0.4]);
set(get(ax_4,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Average duration of food-well visit (s)');
line(repmat(xdata,[2 1]) + control_xoffset, ...
    [min(duration(control_idx,:)); max(duration(control_idx,:))], ...
    'LineStyle',control_style.LineStyle, ...
    'LineWidth',control_style.LineWidth, ...
    'Color',control_style.Color);
line(xdata + control_xoffset, ...
    median(duration(control_idx,:),1),control_style);
line(repmat(xdata,[2 1]) + lesion_xoffset, ...
    [min(duration(lesion_idx,:)); max(duration(lesion_idx,:))], ...
    'LineStyle',lesion_style.LineStyle, ...
    'LineWidth',lesion_style.LineWidth, ...
    'Color',lesion_style.Color);
line(xdata + lesion_xoffset, ...
    median(duration(lesion_idx,:),1),lesion_style);

pause;
print(gcf,'-depsc','-r600','/home/smkim/lesion_manuscript/revisions/linear_track.eps');
delete(gcf);

%}
