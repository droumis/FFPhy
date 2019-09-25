% plots the time taken for each trial over days



Veqn = '>=0'
minV =  str2num(Veqn(end))
maxstage = 3% [1 2 3]
minVPF = 2 %cm/sec
minPeakPF = 3
lessthan=0
includestates = 6

%Animal selection
%-----------------------------------------------------
animals = {'M3'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filter

days='[1:18]';%,'1:10';



%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

epochfilter{1} = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch, 2)'];
%cellfilter = '(isequal($area, ''CA1'') && ($meanrate <7))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '$velocity <0'} };
%timefilter = { {'JY_getriptimes','($nripples ==0)', [], 3,'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { {'JY_getlinvelocity', '(($velocity) >= 0))', 6} };
%f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter);


%only include cells with placefields
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'JY_singleepochanal';

f = setfilteriterator(f,iterator);

out=setfilterfunction(f, 'JY_getreroutetrialtime', {'data'});
        
out=runfilter(out);


% plot trajectory time
% setup variables
winsize=10;
epoch_len = [];
day_len = [];
trialsegments_all = [];
day_totalseg = [];
numseg_mean = [];
numseg_std = [];
day_totaltime = [];
totaltime_mean = [];
totaltime_std = [];
trials_day=[];
barrier=[];


uniquedays=unique(out.epochs{1,1}(:,1));

% loop through days
for i = 1:size(uniquedays)
    trials_day = [];
    jind=find(out.epochs{1,1}(:,1)==uniquedays(i,1));
    % collecting trial segments across epochs and store which trial each day and
    % epoch ends so it can be used to plot later
    for j = 1:size(jind,1);
    trials_day=[trials_day; out.output{1,1}{1,jind(j)}.trialduration];
    % collecting trial segments across epochs and store which trial each day and
    % epoch ends so it can be used to plot later
    day_len = [day_len; out.output{1,1}{1,jind(j)}.ntrial];
    barrier = [barrier; out.output{1,1}{1,jind(j)}.barrier];
   
    end

    epochs = [];
  
    
    % Collect the time of each trial across all epochs in the day
    
    totaltime_mean = [totaltime_mean; moving(trials_day,winsize,@mean)];
    totaltime_std = [totaltime_std; moving(trials_day,winsize,@std)];
    
    %trialsegments_all = [trialsegments_all trialsegments_day];
    
end

epochsummary=[out.epochs{1,1}(:,1) day_len];
totaltime_mean = totaltime_mean/10000; 
totaltime_std = totaltime_std/10000; 


% Plotting time per trial
figure;
set(gcf,'position',[0 0 800 250]);
set(gcf,'PaperPositionMode','auto');
set(gca,'FontSize',12);

% make grey/white strips
o=1;
col=[];

uniquedays=unique(epochsummary(:,1));

while o<size(uniquedays,1)+1;
    currcol=[1 1 1;0.8 0.8 0.8];
    col=[col;currcol];
    o=o+2;
end

daycumsum=[];
for i=1:size(uniquedays,1)
    daycumsum=[daycumsum; sum(epochsummary(epochsummary(:,1)==uniquedays(i,1),2))];
end

daycumsum=[0;cumsum(daycumsum)];

n=0;

for n=1:size(daycumsum,1)-1;
    linen=line([daycumsum(n,1) daycumsum(n+1,1)],[0 0],'Color',col(n,:),'LineWidth',1000);
    text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),90,num2str(uniquedays(n,1)),'HorizontalAlignment','center','FontSize',12);
    set(get(get(linen,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off'); % Exclude line from legend
%     text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),1.15,[num2str(daylength(n+1,2)/100,3) 'm'],'HorizontalAlignment','center','FontSize',10);
    n=n+1;
end
hold on;
% add line to show barrier trials

barrieron=find(barrier(:,1==1));
for i=1:size(barrieron,1);
    lineb=line([barrieron(i) barrieron(i)],[75 75],'LineStyle','s', 'Color','g','LineWidth',2);
    set(get(get(lineb,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off'); % Exclude line from legend
    hold on;
end

plot(totaltime_mean,'-b','LineWidth',2);

plot(totaltime_mean+totaltime_std,':r','LineWidth',2);
plot(totaltime_mean-totaltime_std,':r','LineWidth',2);

% annotating

% for ii = 1:length(day_len)
%     plot([sum(day_len(1:ii)) sum(day_len(1:ii))]-(winsize-1)*(ii), [0 200], 'k');
%     %text(sum(day_len(1:ii))-(winsize)*ii, 85, sprintf('Day %d\nn=%d\ntime=%d',days(ii),day_len(ii),floor(day_totaltime(ii))), ...
%         %'HorizontalAlignment','right','FontSize',10);
%     
% end

xlim([0 max(daycumsum(:,1))]);
ylim([0 100]);
title(sprintf('Trial duration for %s, day %d-%d \n Moving average (%d trial window)', out.animal{1,1},min(uniquedays),max(uniquedays), winsize));
xlabel('Trial');
ylabel('Time (s)');


% % saving

directoryname=f.animal{1,2};
 [s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
 print(sprintf('%sPlot/behav/%s_time',directoryname,f.animal{1,1}),'-depsc');

