
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare dwell time for rewarded versus unrewarded on day 10
%
%{
load('Wtrack_summary');
lesion_idx = find(strcmp({Wtrack_summary(:).group},'lesion'));
control_idx = find(strcmp({Wtrack_summary(:).group},'control'));

for s = 1:length(Wtrack_summary)
    load([Wtrack_summary(s).subject '_Wtrack_journeys']);
    correct_dwelltimes{s} = [ ...
        journeys{2}{1}.visitduration(journeys{2}{1}.correct == 1), ...
        journeys{2}{2}.visitduration(journeys{2}{2}.correct == 1) ];
    correct_dwelltimes{s} = correct_dwelltimes{s}(~isnan(correct_dwelltimes{s}));
    error_dwelltimes{s} = [ ...
        journeys{2}{1}.visitduration(journeys{2}{1}.correct == 0), ...
        journeys{2}{2}.visitduration(journeys{2}{2}.correct == 0) ];
    error_dwelltimes{s} = error_dwelltimes{s}(~isnan(error_dwelltimes{s}));
end

edges = 0:0.5:20;
for s = 1:length(Wtrack_summary)
    n = histc(correct_dwelltimes{s},edges);
    figure(1);
    if strcmp(Wtrack_summary(s).group,'lesion')
        subplot(5,2,s);
        h = bar(edges,n,'histc');
        set(h,'FaceColor','r')
        set(gca,'XLim',[0 20],'YLim',[0 30]);
    elseif strcmp(Wtrack_summary(s).group,'control')
        subplot(5,2,s);
        h = bar(edges,n,'histc');
        set(h,'FaceColor','k')
        set(gca,'XLim',[0 20],'YLim',[0 30]);
    end

    n = histc(error_dwelltimes{s},edges);
    figure(2);
    if strcmp(Wtrack_summary(s).group,'lesion')
        subplot(5,2,s);
        h = bar(edges,n,'histc');
        set(h,'FaceColor','r')
        set(gca,'XLim',[0 20],'YLim',[0 30]);
    elseif strcmp(Wtrack_summary(s).group,'control')
        subplot(5,2,s);
        h = bar(edges,n,'histc');
        set(h,'FaceColor','k')
        set(gca,'XLim',[0 20],'YLim',[0 30]);
    end

end
%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot examples of journeys on W track

trace_style.FaceColor = 'none';
trace_style.EdgeColor = 'flat';
trace_style.Marker = '.';
trace_style.MarkerSize = 10;
trace_style.LineStyle = '-';
trace_style.LineWidth = 3;

axes_style.XLim = [-60 60];
axes_style.YLim = [-20 170];
axes_style.YDir = 'reverse';
axes_style.DataAspectRatio = [1 1 1];
axes_style.Visible = 'off';
axes_style.CLim = [0 120];
colormap('jet');
colorbar('Location','EastOutside');

load('Wtrack_template');
% a control subject
load('M06_Wtrack_smoothedpos');
load('M06_Wtrack_journeys');
% first run on day 10 (fluent performance)
startzone = journeys{10}{1}.startzone;
endzone = journeys{10}{1}.endzone;
correct = journeys{10}{1}.correct;
snippets = journeys{10}{1}.smoothedpos_snippets;

% draw the track
patch('Vertices',Wtrack_template.controlpoints, ...
    'Faces',[1 2 8 9 3 4 10 11 5 6 12 7], ...
    'FaceColor',[1 1 1],'EdgeColor',[0 0 0], ...
    'LineWidth',2);
% draw disc over food-well zone
for z = 1:length(Wtrack_template.zones)
    patch_handle(z) = patch( ...
        'Vertices',Wtrack_template.zones(z).vertices, ...
        'Faces',1:size(Wtrack_template.zones(z).vertices,1), ...
        'FaceColor',[0.8 0.8 0.8],'EdgeColor','none', ...
        'Visible','off');
end
for i = 37:43
    set(patch_handle(endzone(i)),'Visible','on');
    posdata = [ snippets(i).data; nan(1,size(snippets(i).data,2)) ];
    %posdata = snippets(i).data;
    speed = hypot(posdata(:,4),posdata(:,5));
    h = patch('Vertices',posdata(:,2:3), ...
        'Faces',1:size(posdata,1), ...
        'FaceVertexCData',speed, ...
        trace_style);
    set(gca,axes_style);
    print(gcf,'-depsc','-r600', ...
        sprintf('/home/smkim/lesion_manuscript/journey_%02d',i));
    delete(h);
    for z = 1:3
        set(patch_handle(z),'Visible','off');
    end
end
%}

%{
% summary median+range plots across days
load('lineartrack_summary');
load('Wtrack_summary');
days = 1:10;

xdata = [-2.5 -1 NaN 1:10];
axes_style.Color = 'none';
axes_style.XLim = [-3 10.5];
%axes_style.YLim = [0 1];
%axes_style.YTick = 0:0.2:1;
%axes_style.YLim = [0 250];
%axes_style.YTick = 0:50:250;
%axes_style.YLim = [0 90];
%axes_style.YTick = 0:15:90;
axes_style.YLim = [0 18];
axes_style.YTick = 0:3:18;
axes_style.FontSize = 18;
axes_style.FontName = 'Arial';
axes_style.XTick = [-2.5 -1 -0.5 0.5 1:10];
axes_style.LineWidth = 3;
axes_style.TickDir = 'out';
axes_style.TickLength = [0.02 0.02];
axes_style.XTickLabel = ['Pre','Post',' ',' ',arrayfun(@num2str,1:10,'UniformOutput',0)];
axes_style.Position = [0.2 0.2 0.7 0.3];
axes_style.Box = 'on';

lesion_style.Color = 'r';
lesion_style.LineStyle = '-';
lesion_style.Marker = 'o';
lesion_style.LineWidth = 3;
lesion_style.MarkerSize = 8;
lesion_style.MarkerFaceColor = 'w';
control_style.Color = 'k';
control_style.LineStyle = '-';
control_style.Marker = 's';
control_style.LineWidth = 3;
control_style.MarkerSize = 8;
control_style.MarkerFaceColor = 'k';

lesion_idx = strcmp({Wtrack_summary(:).group},'lesion');
control_idx = strcmp({Wtrack_summary(:).group},'control');

lineartrack_pcorr = vertcat(lineartrack_summary(:).correct) ./ ...
    vertcat(lineartrack_summary(:).total);
Wtrack_pcorr = vertcat(Wtrack_summary(:).correct) ./ ...
    vertcat(Wtrack_summary(:).total);

lineartrack_total = vertcat(lineartrack_summary(:).total);
Wtrack_total = vertcat(Wtrack_summary(:).total);

lineartrack_speed = vertcat(lineartrack_summary(:).overallspeed);
Wtrack_speed = vertcat(Wtrack_summary(:).overallspeed);

lineartrack_visitduration = vertcat(lineartrack_summary(:).visitduration);
Wtrack_visitduration = vertcat(Wtrack_summary(:).visitduration);

% median+whisker plots
lesion_median = [ ...
    median(lineartrack_visitduration(lesion_idx,:),1), NaN, ...
    median(Wtrack_visitduration(lesion_idx,:),1) ];
lesion_range = [ ...
    min(lineartrack_visitduration(lesion_idx,:)), NaN, ...
    min(Wtrack_visitduration(lesion_idx,:)); ...
    max(lineartrack_visitduration(lesion_idx,:)), NaN, ...
    max(Wtrack_visitduration(lesion_idx,:)) ];
control_median = [ ...
    median(lineartrack_visitduration(control_idx,:),1), NaN, ...
    median(Wtrack_visitduration(control_idx,:),1) ];
control_range = [ ...
    min(lineartrack_visitduration(control_idx,:)), NaN, ...
    min(Wtrack_visitduration(control_idx,:)); ...
    max(lineartrack_visitduration(control_idx,:)), NaN, ...
    max(Wtrack_visitduration(control_idx,:)) ];

lesion_xoffset = +0.11;
control_xoffset = -0.11;

line(repmat(xdata,[2 1])+control_xoffset,control_range,'LineWidth',3,'Color','k');
control_line = line(xdata+control_xoffset,control_median,control_style);
line(repmat(xdata,[2 1])+lesion_xoffset,lesion_range,'LineWidth',3,'Color','r');
lesion_line = line(xdata+lesion_xoffset,lesion_median,lesion_style);
%{
[legend_h,object_h,plot_h,text_strings] = ...
    legend([control_line lesion_line], ...
    ['Control (n = ' num2str(nnz(control_idx)) ')'], ...
    ['Hippocampal lesion (n = ' num2str(nnz(lesion_idx)) ')'], ...
    'Location','BestOutside');
set(object_h(1),'FontSize',18,'FontName','Arial');
set(object_h(2),'FontSize',18,'FontName','Arial');
legend(gca,'boxoff');
line([-3 -3 -0.5 -0.5],[0.4 0.45 0.45 0.4],'LineWidth',3,'Color','k');
line([0.5 0.5 10.5 10.5],[0.4 0.45 0.45 0.4],'LineWidth',3,'Color','k');
%}

%set(get(gca,'YLabel'),'FontSize',18,'FontName','Arial','Color','k','String','Proportion correct')
%set(get(gca,'YLabel'),'FontSize',18,'Color','k','String','Number of journeys');
%set(get(gca,'YLabel'),'FontSize',18,'Color','k','String','Running speed')
set(get(gca,'YLabel'),'FontSize',18,'Color','k','String','Visit duration')

%set(get(gca,'XLabel'),'FontSize',18,'FontName','Arial','Color','k','String','Day');
set(gca,axes_style);
%}

%{
%%%%% plot deficit in visiting center arm on day 1
%
axes_style.XLim = [0 0.55];
axes_style.XTick = 0:0.1:0.5;
axes_style.YLim = [0.63 1];
axes_style.YTick = [0.67, 0.7:0.1:1];
axes_style.FontSize = 18;
axes_style.FontName = 'Arial';
axes_style.LineWidth = 3;
axes_style.TickDir = 'out';
axes_style.TickLength = [0.03 0.03];
axes_style.Box = 'off';
axes_style.Layer = 'bottom';
axes_style.Color = 'none';
axes_style.DataAspectRatio = [1 1 1];
axes_style.Visible = 'on';
axes_style.YDir = 'normal';

lesion_style.Color = 'r';
lesion_style.LineStyle = 'none';
lesion_style.Marker = 'o';
lesion_style.LineWidth = 3;
lesion_style.MarkerSize = 8;
lesion_style.MarkerFaceColor = 'none';
control_style.Color = 'k';
control_style.LineStyle = 'none';
control_style.Marker = 's';
control_style.LineWidth = 3;
control_style.MarkerSize = 8;
control_style.MarkerFaceColor = 'k';

load('Wtrack_summary');
lesion_idx = find(strcmp({Wtrack_summary(:).group},'lesion'));
control_idx = find(strcmp({Wtrack_summary(:).group},'control'));
p_center = vertcat(Wtrack_summary(:).center) ./ ...
    vertcat(Wtrack_summary(:).total);
p_choicepoint = vertcat(Wtrack_summary(:).choicepoint) ./ ...
    vertcat(Wtrack_summary(:).total);

for day = 1:2;
    figure();
    axes('Position',[0.1 0.1 0.55 0.5],axes_style);
    % plot guidelines for optimal task performance
    line(get(gca,'XLim'),[1 1],'LineStyle','--','LineWidth',1,'Color','k');
    line([0.5 0.5],get(gca,'YLim'),'LineStyle','--','LineWidth',1,'Color','k');
    % plot scatterplot
    line(p_center(control_idx,day),p_choicepoint(control_idx,day),control_style);
    line(p_center(lesion_idx,day),p_choicepoint(lesion_idx,day),lesion_style);
    %print(gcf,'-depsc','-r600',sprintf('/home/smkim/lesion_manuscript/journey%d',day));
    delete(gcf);
end
%}

%{
% plot examples of side-to-side journeys on W track
figure(2);
trace_style.FaceColor = 'none';
trace_style.EdgeColor = 'flat';
trace_style.Marker = '.';
trace_style.MarkerSize = 10;
trace_style.LineStyle = '-';
trace_style.LineWidth = 3;

axes_style.XLim = [-60 60];
axes_style.YLim = [-20 170];
axes_style.YDir = 'reverse';
axes_style.DataAspectRatio = [1 1 1];
axes_style.Visible = 'off';
axes_style.CLim = [0 120];
colormap('jet');
colorbar('Location','EastOutside');

load('Wtrack_template');
% a lesion subject
load('M12_Wtrack_journeys');
% first run on day 10 (fluent performance)
startzone = journeys{1}{1}.startzone;
endzone = journeys{1}{1}.endzone;
correct = journeys{1}{1}.correct;
snippets = journeys{1}{1}.smoothedpos_snippets;

% draw the track
patch('Vertices',Wtrack_template.controlpoints, ...
    'Faces',[1 2 8 9 3 4 10 11 5 6 12 7], ...
    'FaceColor',[1 1 1],'EdgeColor',[0 0 0], ...
    'LineWidth',2);
for i = 19:23
    posdata = [ snippets(i).data; nan(1,size(snippets(i).data,2)) ];
    %posdata = snippets(i).data;
    speed = hypot(posdata(:,4),posdata(:,5));
    h = patch('Vertices',posdata(:,2:3), ...
        'Faces',1:size(posdata,1), ...
        'FaceVertexCData',speed, ...
        trace_style);
    set(gca,axes_style);
    print(gcf,'-depsc','-r600',sprintf('/home/smkim/lesion_manuscript/journey_%02d',i));
    delete(h);
end
%}


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot example learning curves
clear('axes_style');
load('learning_curves_2');
lesion_idx = find(strcmp({learning_curves(:).group},'lesion'));
control_idx = find(strcmp({learning_curves(:).group},'control'));

%for s = 1:2
for s = 1:length(learning_curves)
    x = learning_curves(s);

    % day 1 consists of journeys daybounds(1) to daybounds(2)-1, inclusive;
    % day 2 consists of journeys daybounds(2) to daybounds(3)-1, inclusive; etc.
    daybounds = [1; find(diff(x.dayrun(:,1))); size(x.dayrun,1)];

    % axes boundaries
    max_count = length(x.prob);
    criterion = x.criterion;
    axes_style.XLim = [-10 250*ceil(max_count/250)];
    axes_style.YLim = [0 1];
    axes_style.FontSize = 14;
    axes_style.FontName = 'Arial';
    axes_style.XTick = [1 250:250:axes_style.XLim(2)];
    axes_style.YTick = [0 0.25 0.5 0.75 1];
    axes_style.LineWidth = 2;
    axes_style.TickDir = 'out';
    axes_style.TickLength = [0.02 0.02]*2010/diff(axes_style.XLim);
    axes_style.Box = 'on';
    axes_style.Layer = 'top';
    axes_style.Color = 'none';

    point_style.Color = 'k';
    point_style.LineStyle = 'none';
    point_style.Marker = 'o';
    point_style.MarkerFaceColor = 'k';
    point_style.MarkerSize = 1;

    whisker_style.Color = [0.75 0.75 0.75];
    whisker_style.LineStyle = '-';
    whisker_style.LineWidth = 1;
    whisker_style.Marker = 'none';

    % draw reference patches for demarcating boundaries between days
    a3 = axes('Position',[0.2 0.075 0.7*diff(axes_style.XLim)/2010 0.85],axes_style);
    set(a3,'Visible','off');
    axes(a3);
    for i = 2:length(daybounds)
        if mod(i,2)
            lightgreen = hex2dec({'e0' 'fa' 'ff'})'/255;
            patch([daybounds(i-1) daybounds(i-1) daybounds(i) daybounds(i)], ...
                [0 1 1 0],lightgreen);
            text(0.5*(daybounds(i-1)+daybounds(i)),0.5,num2str(i-1), ...
                'Color','k','FontSize',14,'FontName','Arial', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        else
            lightblue = hex2dec({'ff' 'ff' 'e0'})'/255;
            patch([daybounds(i-1) daybounds(i-1) daybounds(i) daybounds(i)], ...
                [0 1 1 0],lightblue);
            text(0.5*(daybounds(i-1)+daybounds(i)),0.5,num2str(i-1), ...
                'Color','k','FontSize',14,'FontName','Arial', ...
                'HorizontalAlignment','center','VerticalAlignment','middle');
        end
    end
    
    a1 = axes('Position',[0.2 0.55 0.7*diff(axes_style.XLim)/2010 0.35],axes_style);
    set(a1,'XAxisLocation','top','YLim',[0 1]);
    set(get(a1,'YLabel'),'FontSize',14,'Color','k','String','Inbound');
    set(get(a1,'XLabel'),'String','Cumulative count of journeys', ...
        'FontSize',14,'FontName','Arial','Color','k');

    a2 = axes('Position',[0.2 0.1 0.7*diff(axes_style.XLim)/2010 0.35],axes_style);
    set(a2,'XAxisLocation','bottom','YLim',[0 1]);
    set(get(a2,'YLabel'),'FontSize',14,'Color','k','String','Outbound');
    set(get(a2,'XLabel'),'String','Cumulative count of journeys', ...
        'FontSize',14,'FontName','Arial','Color','k');

    % draw whiskers
    for j = 1:numel(x.correct)
        if strcmp(x.type{j},'outbound')
            axes(a2);
            line([j j],x.prob(j,2:3),whisker_style);
        elseif strcmp(x.type{j},'inbound')
            axes(a1);
            line([j j],x.prob(j,2:3),whisker_style);
        end
    end
    % draw points for outbound
    outpick = find(strcmp(x.type,'outbound'));
    axes(a2);
    line(outpick,x.prob(outpick,1),point_style);
    % draw line for inbound
    inpick = find(strcmp(x.type,'inbound'));
    axes(a1);
    line(inpick,x.prob(inpick,1),point_style);

    % redraw point+whisker in red on learning trial, if it exists
    if ~isnan(x.lt_outbound)
        axes(a2);
        ld_out = x.ld_outbound;
        lt_out = x.lt_outbound;
        line([lt_out lt_out], ...
            x.prob(lt_out,2:3),whisker_style,'Color','r');
        line(lt_out,x.prob(lt_out,1),point_style, ...
            'Color',[0.75 0 0],'MarkerFaceColor',[0.75 0 0],'MarkerSize',3);
        text(lt_out,x.prob(lt_out,2)-0.02, ...
            ['day ' num2str(ld_out)], ...
            'FontSize',14,'FontName','Arial','Color','r','VerticalAlignment','top');
    end
    if ~isnan(x.lt_inbound)
        axes(a1);
        ld_in = x.ld_inbound;
        lt_in = x.lt_inbound;
        line([lt_in lt_in], ...
            x.prob(lt_in,2:3),whisker_style,'Color','r');
        line(lt_in,x.prob(lt_in,1),point_style, ...
            'Color',[0.75 0 0],'MarkerFaceColor',[0.75 0 0],'MarkerSize',3);
        text(lt_in,x.prob(lt_in,2)-0.02, ...
            ['day ' num2str(ld_in)], ...
            'FontSize',14,'FontName','Arial','Color','r','VerticalAlignment','top');
    end

    % dashed lines to indicate background probability
    axes(a1);
    line(get(gca,'XLim'),[criterion criterion],'Color','k','LineWidth',1,'LineStyle','--');
    axes(a2);
    line(get(gca,'XLim'),[criterion criterion],'Color','k','LineWidth',1,'LineStyle','--');

    print(gcf,'-depsc','-r600',['/home/smkim/lesion_manuscript/' x.subject '_' x.group '_half.eps']);
    pause;
    delete(gcf);
end
%}

%{
%%%%%%%%%%%%%%
load('learning_curves_2');
clear('axes_style');
% plot days to reach learning criterion on inbound and outbound components
lesion_idx = find(strcmp({learning_curves(:).group},'lesion'));
control_idx = find(strcmp({learning_curves(:).group},'control'));
axes_style.YLim = [0.5 2.5];
axes_style.XLim = [0.5 13];
axes_style.XTick = [1:10, 10.5, 11.5, 12.5];
axes_style.YTick = 1:2;
axes_style.YTickLabel = {'Control','Hippocampal lesion'};
axes_style.FontSize = 18;
axes_style.FontName = 'Arial';
axes_style.LineWidth = 3;
axes_style.TickDir = 'out';
axes_style.TickLength = [0.025 0.025];
axes_style.Box = 'off';
axes_style.Layer = 'top';
axes_style.Color = 'none';
axes_style.YDir = 'reverse';

lesion_style.Color = 'r';
lesion_style.LineStyle = 'none';
lesion_style.Marker = 'o';
lesion_style.LineWidth = 3;
lesion_style.MarkerSize = 8;
lesion_style.MarkerFaceColor = 'w';
control_style.Color = 'k';
control_style.LineStyle = 'none';
control_style.Marker = 's';
control_style.LineWidth = 3;
control_style.MarkerSize = 8;
control_style.MarkerFaceColor = 'k';

ld_out = [learning_curves(:).ld_outbound];
ld_out(isnan(ld_out)) = 12.5;
ld_in = [learning_curves(:).ld_inbound];
ld_in(isnan(ld_in)) = 12.5;
ld_in_ydata(control_idx) = 1*ones(numel(control_idx),1);
ld_in_ydata(lesion_idx) = 2*ones(numel(lesion_idx),1);
ld_out_ydata(control_idx) = 1*ones(numel(control_idx),1);
ld_out_ydata(lesion_idx) = 2*ones(numel(lesion_idx),1);
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

a1 = axes('Position',[0.2 0.6 0.7 0.2],axes_style);
line(ld_in(control_idx),ld_in_ydata(control_idx),control_style);
line(ld_in(lesion_idx),ld_in_ydata(lesion_idx),lesion_style);
a2 = axes('Position',[0.2 0.1 0.7 0.2],axes_style);
line(ld_out(control_idx),ld_out_ydata(control_idx),control_style);
line(ld_out(lesion_idx),ld_out_ydata(lesion_idx),lesion_style);

% ranksum test on days to reach learning criterion
disp(sprintf('ranksum test:   inbound p %f   outbound p %f', ...
    [ ranksum(ld_in(lesion_idx),ld_in(control_idx)) ...
    ranksum(ld_out(lesion_idx),ld_out(control_idx)) ]));
% BF rank test (Brunner)
%in_stats = BFranktest(ld_in(lesion_idx),ld_in(control_idx),'permutation',1e5);
%out_stats = BFranktest(ld_out(lesion_idx),ld_out(control_idx),'permutation',1e5);
%disp(sprintf('BF rank test:   inbound p %f   outbound p %f',[in_stats.p, out_stats.p]));
%
%}

%{
%%%%%%%%%%%%%%
% plot trials to reach learning criterion on inbound and outbound components
load('behavperform_3');
clear('axes_style');
lesion_idx_learned_inbound = find(strcmp({behavperform(:).group},'lesion') & ...
    ([behavperform(:).lt_in] > 0));
lesion_idx_failed_inbound = find(strcmp({behavperform(:).group},'lesion') & ...
    ([behavperform(:).lt_in] < 0));
lesion_idx_learned_outbound = find(strcmp({behavperform(:).group},'lesion') & ...
    ([behavperform(:).lt_out] > 0));
lesion_idx_failed_outbound = find(strcmp({behavperform(:).group},'lesion') & ...
    ([behavperform(:).lt_out] < 0));

control_idx = find(strcmp({behavperform(:).group},'control'));
axes_style.XLim = [0 750];
axes_style.XTick = 0:150:750;
axes_style.YLim = [0.5 2.5];
axes_style.YTick = 1:2;
axes_style.YTickLabel = {'Control','Hippocampal lesion'};
axes_style.FontSize = 18;
axes_style.FontName = 'Arial';
axes_style.LineWidth = 3;
axes_style.TickDir = 'out';
axes_style.TickLength = [0.025 0.025];
axes_style.Box = 'off';
axes_style.Layer = 'top';
axes_style.Color = 'none';
axes_style.YDir = 'reverse';

lesion_style.Color = 'r';
lesion_style.LineStyle = 'none';
lesion_style.MarkerSize = 8;
lesion_style.MarkerFaceColor = 'w';
control_style.Color = 'k';
control_style.LineStyle = 'none';
control_style.Marker = 's';
control_style.MarkerSize = 8;
control_style.MarkerFaceColor = 'k';

lt_in = abs([behavperform(:).lt_in]);
lt_in_ydata(control_idx) = 1*ones(numel(control_idx),1) ;
lt_in_ydata(lesion_idx) = 2*ones(numel(lesion_idx),1);

lt_out = abs([behavperform(:).lt_out]);
lt_out_ydata(control_idx) = 1*ones(numel(control_idx),1);
lt_out_ydata(lesion_idx) = 2*ones(numel(lesion_idx),1);

% jitter points with negative-valued learning trials (which indicate that they
% did not learn) 
jitter = 0.1;
lt_in_ydata = lt_in_ydata + jitter*([behavperform(:).lt_in] < 0);
lt_out_ydata = lt_out_ydata + jitter*([behavperform(:).lt_out] < 0);

a1 = axes('Position',[0.2 0.6 0.7 0.2],axes_style);
line(lt_in(control_idx),lt_in_ydata(control_idx),control_style);
line(lt_in(lesion_idx_learned_inbound), ...
    lt_in_ydata(lesion_idx_learned_inbound),lesion_style,'Marker','o');
line(lt_in(lesion_idx_failed_inbound), ...
    lt_in_ydata(lesion_idx_failed_inbound),lesion_style,'Marker','x');
xlabel('Number of inbound trials to reach criterion');
a2 = axes('Position',[0.2 0.1 0.7 0.2],axes_style);
line(lt_out(control_idx),lt_out_ydata(control_idx),control_style);
line(lt_out(lesion_idx_learned_outbound), ...
    lt_out_ydata(lesion_idx_learned_outbound),lesion_style,'Marker','o');
line(lt_out(lesion_idx_failed_outbound), ...
    lt_out_ydata(lesion_idx_failed_outbound),lesion_style,'Marker','x');
xlabel('Number of outbound trials to reach criterion');

% ranksum test on trials to reach learning criterion
disp(sprintf('ranksum test:   inbound p %f   outbound p %f', ...
    [ ranksum(lt_in(lesion_idx),lt_in(control_idx)) ...
    ranksum(lt_out(lesion_idx),lt_out(control_idx)) ]));
% BF rank test (Brunner)
%in_stats = BFranktest(lt_in(lesion_idx),lt_in(control_idx),'permutation',1e5);
%out_stats = BFranktest(lt_out(lesion_idx),lt_out(control_idx),'permutation',1e5);
%disp(sprintf('BF rank test:   inbound p %f   outbound p %f',[in_stats.p, out_stats.p]));
%
%}

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
      0.85*diff(in_axes_style.XLim)/(max_trials_anywhere+left_margin) ...
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
      in_axes_style.Position(3) ...
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

  %print(gcf,'-depsc','-r600','/home/smkim/lesion_manuscript/revisions/test.eps');
  print(gcf,'-depsc','-r600',sprintf('/home/smkim/lesion_manuscript/revisions/%s_criterion(%f).eps',subject_label,x.criterion));
  pause;
  delete(gcf);
end
%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot example MA curves, separately for inbound and outbound
clear('axes_style');
load('behavperform_3');
lesion_idx = find(strcmp({behavperform(:).group},'lesion'));
control_idx = find(strcmp({behavperform(:).group},'control'));
% maximum number of trials on either inbound or outbound for any subject
max_trials_anywhere = 1100;

MA_order = 40;

FONTSIZE = 14;

for s = 1:length(behavperform)
  x = behavperform(s);
  criterion = x.criterion;
  max_in = length(x.inreward);
  max_out = length(x.outreward);
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
      0.85*diff(in_axes_style.XLim)/(max_trials_anywhere+left_margin) ...
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
      in_axes_style.Position(3) ...
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

  %print(gcf,'-depsc','-r600','/home/smkim/lesion_manuscript/revisions/test.eps');
  print(gcf,'-depsc','-r600',sprintf('/home/smkim/lesion_manuscript/revisions/%s_MA%d.eps',subject_label,MA_order));
  pause;
  delete(gcf);
end
%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tabulate various things about the learning curves, 
% separately for inbound and outbound, for every subject, 
% with between-group p-values:
%   total number of trials
%   trials to reach criterion
%   day on which learning criterion was achieved
%   slope of learning curve at criterion
%   estimated probability correct at the end of testing
%   
%   
%   
clear;
load('behavperform_2');
lesion_idx = find(strcmp({behavperform(:).group},'lesion'));
control_idx = find(strcmp({behavperform(:).group},'control'));

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
  learning_stats(s).final_in = x.inprobcorrect(end,1);
  if ~isnan(x.ld_in)
    % select window around the learning trial
    window_in = ((x.lt_in-halfwidth):(x.lt_in+halfwidth))';
    % approximate inverse variances; for a Gaussian distribution, 95% confidence
    % intervals covers the middle 4-sigma of the data
    weights = ((x.inprobcorrect(window_in,3) - ...
        x.inprobcorrect(window_in,2))/4) .^ (-2);
    coeffs_in = lscov([window_in, ones(size(window_in))], ...
        x.inprobcorrect(window_in,1),weights);
    learning_stats(s).slope_in = coeffs_in(1);
    learning_stats(s).lt_in = x.lt_in;
  else
    learning_stats(s).slope_in = NaN;
    learning_stats(s).lt_in = NaN;
  end

  learning_stats(s).trialcount_out = length(x.outreward);
  learning_stats(s).ld_out = x.ld_out;
  learning_stats(s).final_out = x.outprobcorrect(end,1);
  if ~isnan(x.ld_out)
    % select window around the learning trial
    window_out = ((x.lt_out-halfwidth):(x.lt_out+halfwidth))';
    % approximate inverse variances; for a Gaussian distribution, 95% confidence
    % intervals covers the middle 4-sigma of the data
    weights = ((x.outprobcorrect(window_out,3) - ...
        x.outprobcorrect(window_out,2))/4) .^ (-2);
    coeffs_out = lscov([window_out, ones(size(window_out))], ...
        x.outprobcorrect(window_out,1),weights);
    learning_stats(s).slope_out = coeffs_out(1);
    learning_stats(s).lt_out = x.lt_out;
  else
    learning_stats(s).slope_out = NaN;
    learning_stats(s).lt_out = NaN;
  end
end

FONTSIZE = 14;

lesion_style.Color = 'r';
lesion_style.LineStyle = 'none';
lesion_style.Marker = 'o';
lesion_style.LineWidth = 3;
lesion_style.MarkerSize = 8;
lesion_style.MarkerFaceColor = 'w';
control_style.Color = 'k';
control_style.LineStyle = 'none';
control_style.Marker = 's';
control_style.LineWidth = 3;
control_style.MarkerSize = 8;
control_style.MarkerFaceColor = 'k';

% plot days to reach learning criterion on inbound and outbound components
axes_style.YLim = [0.5 2.5];
axes_style.XLim = [0.5 13];
axes_style.XTick = [1:10, 10.5, 11.5, 12.5];
axes_style.XTickLabel = horzcat(arrayfun(@num2str,1:10, ...
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
ld_out(isnan(ld_out)) = 12.5;
ld_in = [learning_stats(:).ld_in];
ld_in(isnan(ld_in)) = 12.5;
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
a1 = axes('Position',[0.3 0.7 0.5 0.2],axes_style);
set(get(a1,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Days to reach learning criterion');
text(6.5,2.5,'Inbound','FontSize',FONTSIZE,'FontName','Arial', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center');
line(ld_in(control_idx),ld_in_ydata(control_idx),control_style);
line(ld_in(lesion_idx),ld_in_ydata(lesion_idx),lesion_style);
a2 = axes('Position',[0.3 0.3 0.5 0.2],axes_style);
set(get(a2,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Days to reach learning criterion');
text(6.5,2.5,'Outbound','FontSize',FONTSIZE,'FontName','Arial', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center');
line(ld_out(control_idx),ld_out_ydata(control_idx),control_style);
line(ld_out(lesion_idx),ld_out_ydata(lesion_idx),lesion_style);
% ranksum test on days to reach learning criterion
disp(sprintf('days to criterion, ranksum test:   inbound p %f   outbound p %f', ...
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

final_out = [learning_stats(:).final_out];
final_in = [learning_stats(:).final_in];
final_in_ydata(control_idx) = 2*ones(numel(control_idx),1);
final_in_ydata(lesion_idx) = 1*ones(numel(lesion_idx),1);
final_out_ydata(control_idx) = 2*ones(numel(control_idx),1);
final_out_ydata(lesion_idx) = 1*ones(numel(lesion_idx),1);
jitter = 0.15;
final_out_ydata = final_out_ydata + jitter*randn(size(final_out_ydata));
final_in_ydata = final_in_ydata + jitter*randn(size(final_in_ydata));
figure(2);
a1 = axes('Position',[0.3 0.7 0.5 0.2],axes_style);
set(get(a1,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Final estimated probability correct');
text(0.5,2.5,'Inbound','FontSize',FONTSIZE,'FontName','Arial', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center');
line(final_in(control_idx),final_in_ydata(control_idx),control_style);
line(final_in(lesion_idx),final_in_ydata(lesion_idx),lesion_style);
a2 = axes('Position',[0.3 0.3 0.5 0.2],axes_style);
set(get(a2,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Final estimated probability correct');
text(0.5,2.5,'Outbound','FontSize',FONTSIZE,'FontName','Arial', ...
    'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center');
line(final_out(control_idx),final_out_ydata(control_idx),control_style);
line(final_out(lesion_idx),final_out_ydata(lesion_idx),lesion_style);
% ranksum test on asymptotic performance
disp(sprintf('asymptotic performance, ranksum test:   inbound p %f   outbound p %f', ...
    [ ranksum(final_in(lesion_idx),final_in(control_idx)) ...
    ranksum(final_out(lesion_idx),final_out(control_idx)) ]));
print(gcf,'-depsc','-r600','/home/smkim/lesion_manuscript/revisions/final_performance.eps');
%
pause;
delete(1); delete(2);
%}

%
%
% analysis of errors on inbound trials
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
  inbound_summary(s) = struct(...
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
            outbound_trials(end+1,1) = struct( ...
              'origin',journeys{d}{r}.startzone(j), ...
              'destination',journeys{d}{r}.endzone(j), ...
              'correct',journeys{d}{r}.correct(j), ...
              'running_speed',journeys{d}{r}.speed(j), ...
              'visit_duration',journeys{d}{r}.visitduration(j) );
          else
            inbound_trials(end+1,1) = struct( ...
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
%

% lower-triangle plot, with return errors versus side-to-side errors on the two
% axes
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
lesion_style.LineWidth = 3;
lesion_style.MarkerSize = 8;
lesion_style.MarkerFaceColor = 'none';
control_style.Color = 'k';
control_style.LineStyle = 'none';
control_style.Marker = 's';
control_style.LineWidth = 3;
control_style.MarkerSize = 8;
control_style.MarkerFaceColor = 'k';

a1 = axes(axes_style,'Position',[0.1 0.4 0.35 0.35]);
set(get(a1,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Proportion side-to-side errors');
set(get(a1,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Proportion turn-around errors');
day = 1;
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

a2 = axes(axes_style,'Position',[0.55 0.4 0.35 0.35]);
set(get(a2,'XLabel'),'FontSize',FONTSIZE,'FontName','Arial', ...
    'Color','k','String','Proportion side-to-side errors');
set(get(a2,'YLabel'),'FontSize',FONTSIZE,'FontName','Arial','Color','k', ...
    'String','Proportion turn-around errors');
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

%

% TODO: plot inbound and outbound Wtrack performance separately
