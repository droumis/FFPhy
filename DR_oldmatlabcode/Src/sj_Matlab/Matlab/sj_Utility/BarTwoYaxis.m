function [figHandle] = BarTwoYaxis(bar1,bar2,width,groupLabels,xlabel,y1label,y2label,errors)
% [figHandle] =
% BarTwoYaxis(bar1,bar2,width,groupLabels,xlabel,y1label,y2label,errors)
% plots a 2 axis bar plot by normalizing values to be from 0 to 1 on both
% data sets, plotting them as one grouped bargraph, adding 2nd yaxis
% and then modfying y axis ticks. bar1 & bar2 are column vectors.
% If errors is not empty (Expected 2 row vector same size as numbars) than
% plot errorbars using errorbar function
%
% Author: Shai Shen-Orr
% Date: April 22nd, 2006

mb1 = max(bar1);
mb2 = max(bar2);
nbar1 = bar1/mb1;
nbar2 = bar2/mb2;
numgroups = 2; % number of groups
numbars = length(bar1);
h = bar([nbar1 nbar2],'grouped');

ax1 = gca;
set(ax1,'XColor','k','YColor','b')
set(get(ax1,'Ylabel'),'String',y1label);
set(get(ax1,'Ylabel'),'Color','k');

ax2 = axes('Position',get(ax1,'Position'),...
    'Color','none',...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'XColor','k','YColor','r');
set(ax2,'XTick',[]);
set(get(ax2,'Ylabel'),'String',y2label);
set(get(ax2,'Ylabel'),'Color','k');

ax = get(gcf,'Children');
set(gcf,'CurrentAxes',ax(1));
if(isempty(groupLabels)) 
    groupLabels = [1:numGroups];
end;
axis([0 length(groupLabels)+1 0 1.1]);  
set(gcf,'CurrentAxes',ax(2));
axis([0 length(groupLabels)+1 0 1.1]);  

groupwidth = 0.8;

if(~isempty(errors))
    barvalues = [nbar1 nbar2]';
    xS1 =  get(get(h(1),'Children'),'XData');
    xS2 =  get(get(h(2),'Children'),'XData');
    nerrors = errors;
    nerrors(1,:) = nerrors(1,:)/mb1;
    nerrors(2,:) = nerrors(2,:)/mb2;
    set(gcf,'CurrentAxes',ax1)
    hold on;
    for i = 1:numbars
            x1 = mean(xS1([1 3],i));
            x2 = mean(xS2([1 3],i));
            errorbar([x1 x2], barvalues(:,i), nerrors(:,i), 'k', 'linestyle', 'none');
    end
end;
%Change tick labels
ax1ymax = mb1*1.1;
set(ax1,'YTick',get(ax1,'YTick')*1.1);
set(ax1,'YTickLabel',[0:ax1ymax/(length(get(ax1,'YTick'))-1):ax1ymax]);
ax2ymax = mb2*1.1;
set(ax2,'YTick',get(ax2,'YTick')*1.1);
set(ax2,'YTickLabel',[0:ax2ymax/(length(get(ax2,'YTick'))-1):ax2ymax]);
if(~isempty(groupLabels))
    set(ax1,'XTickLabel',groupLabels);
end;
figHandle = gcf;


