% plots the combined normalised place fields of all places cells active in
% a rippple



Veqn = '>=0'
minV =  str2num(Veqn(end))
maxstage = 3% [1 2 3]
minVPF = 2 %cm/sec
minPeakPF = 3
lessthan=0
includestates = 6

%Animal selection
%-----------------------------------------------------
%animals = {'M3'};

animals = {'K3','L3','M2', 'N3', 'N1','L2','M3','M1','P2'};

%animals = {'M2','P2'};
%animals = {'M2','M1','M3','K3','L2','L3','N1','P2','N3'};

%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filter

days='[1:8]';%,'1:10';

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

outall=setfilterfunction(f, 'JY_gettrajectorydistance', {'data','linpos'});
        
outall=runfilter(outall);

daydistance_summary=[];
for mm=1:size(outall,2)
    
out=outall(1,mm);
uniquedays=unique(out.epochs{1,1}(:,1));
outdata=out.output{1,1};


% plot trajectory distance


% col=jet(length(out.output{1,1}));
% 
% for i=1:length(out.output{1,1});
%     plot(out.output{1,1}{1,i}.trajectorydistance,'.','Color',col(i,:));
%     hold on;
% end
% 
% 
% figure;
% 
% for i=1:length(out.output{1,1});
%     X=[1:1:length(out.output{1,1}{1,i}.normdist)]';
%     Y=out.output{1,1}{1,i}.normdist;
%     Z=ones(length(out.output{1,1}{1,i}.normdist),1)*(length(out.output{1,1})-i+1);
%     plot3(X,Z,Y,'.','Color',col(i,:));
%     %bar3(Y,,'detached');
%     %display(max(out.output{1,1}{1,i}.normdist))
%     grid on;
%     hold on;
% end


% plot trajectory time
% setup variables
winsize=11;
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
daydistance=[];
% loop though data for each day


for i = 1:size(uniquedays)
    trials_day = [];
    jind=find(out.epochs{1,1}(:,1)==uniquedays(i,1));
    
    daydistance(i)=0;
    
    for j = 1:size(jind,1);
    trials_day=[trials_day; out.output{1,1}{1,jind(j)}.normdist];
    % collecting trial segments across epochs and store which trial each day and
    % epoch ends so it can be used to plot later
    day_len = [day_len; out.output{1,1}{1,jind(j)}.ntrial];
    barrier = [barrier; out.output{1,1}{1,jind(j)}.barrier];
    daydistance(i)=daydistance(i)+sum(out.output{1,1}{1,jind(j)}.trajectorydisplacement);
    %daydistance(i)=daydistance(i)+sum(out.output{1,1}{1,jind(j)}.trajectorydistance);
    end
    totaltime_mean =[totaltime_mean; moving(trials_day,winsize,@mean)];
    totaltime_std = [totaltime_std; moving(trials_day,winsize,@std)];
    
    epochs = [];
    
    
    % Collect the time of each trial across all epochs in the day
    
%     totaltime_mean = [totaltime_mean movingstat(trial_time',winsize,@mean)'];
%     totaltime_std = [totaltime_std movingstat(trial_time',winsize,@std)'];
    
    %trialsegments_all = [trialsegments_all trialsegments_day];
    
end

 
epochsummary=[out.epochs{1,1}(:,1) day_len];
 
% % Plotting time per trial
% figure;
% set(gcf,'position',[0 0 800 400]);
% set(gcf,'PaperPositionMode','auto');
% set(gca,'FontSize',12);
% 
% % make grey/white strips
% o=1;
% col=[];
% 
% while o<size(uniquedays,1)+1;
%     currcol=[1 1 1;0.8 0.8 0.8];
%     col=[col;currcol];
%     o=o+2;
% end
% 
% daycumsum=[];
% for i=1:size(uniquedays,1)
%     daycumsum=[daycumsum; sum(epochsummary(epochsummary(:,1)==uniquedays(i,1),2))];
% end
% 
% daycumsum=[0;cumsum(daycumsum)];
% 
% n=0;
% 
% for n=1:size(daycumsum,1)-1;
%     linen=line([daycumsum(n,1) daycumsum(n+1,1)],[0 0],'Color',col(n,:),'LineWidth',1000);
%     text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),6,num2str(uniquedays(n,1)),'HorizontalAlignment','center','FontSize',12);
%     set(get(get(linen,'Annotation'),'LegendInformation'),...
%         'IconDisplayStyle','off'); % Exclude line from legend
% %     text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),1.15,[num2str(daylength(n+1,2)/100,3) 'm'],'HorizontalAlignment','center','FontSize',10);
%     n=n+1;
% end
% hold on;
% %line([1 size(barrier,1)],[1 1],'LineStyle','-', 'Color','k','LineWidth',2);
% 
% % add line to show barrier trials
% 
% barrieron=find(barrier(:,1==1));
% for i=1:size(barrieron,1);
%     lineb=line([barrieron(i) barrieron(i)],[5 5],'LineStyle','s', 'Color','g','LineWidth',2);
%     set(get(get(lineb,'Annotation'),'LegendInformation'),...
%         'IconDisplayStyle','off'); % Exclude line from legend
%     hold on;
% end
% 
% plot(totaltime_mean,'LineWidth',2);
% 
% plot(totaltime_mean+totaltime_std,'r:','LineWidth',2);
% plot(totaltime_mean-totaltime_std,'r:','LineWidth',2);
% 
% % annotating
% 
% % for ii = 1:length(day_len)
% %     plot([sum(day_len(1:ii)) sum(day_len(1:ii))]-(winsize-1)*(ii), [0 200], 'k');
% %     %text(sum(day_len(1:ii))-(winsize)*ii, 85, sprintf('Day %d\nn=%d\ntime=%d',days(ii),day_len(ii),floor(day_totaltime(ii))), ...
% %         %'HorizontalAlignment','right','FontSize',10);
% %     
% % end
% 
% xlim([0 max(daycumsum(:,1))]);
% ylim([0 5]); %max(totaltime_mean+totaltime_std)]);
% title(sprintf('Trial distance for %s, day %d-%d \n Moving average (%d trial window)', out.animal{1,1},min(uniquedays),max(uniquedays), winsize));
% xlabel('Trial');
% ylabel('Multiples of shortest distance');
% 
% [s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
% print(sprintf('%sPlot/behav/%s_distance',directoryname,f.animal{1,1}),'-depsc');


% figure;
% 
% bar(daydistance./100);
% 
% title(sprintf('Distance travelled for %s, day %d-%d', out.animal{1,1},min(uniquedays),max(uniquedays)));
% xlabel('Day');
% ylabel('Distance (m)');

%clear all;

% % saving
% 
% directoryname=f.animal{1,2};
% [s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
% print(sprintf('%sPlot/behav/%s_totaldistance',directoryname,f.animal{1,1}),'-depsc');

daydistance_summary=[daydistance_summary; daydistance];
end


% convert from cm to m
 daydistance_summary=daydistance_summary/100;
 
% 
% % plot raw data
% 
% figure;
% 
% circlesize=100;
% 
% %daydistance_summary(4,:)=daydistance_summary(4,:)/2;
% 
% for mm=1:4
% 
% scatter(1:8,daydistance_summary(mm,:),circlesize,'b','MarkerEdgeColor','b','LineWidth',3); hold on;
% 
% end
% 
% for mm=5:9
% 
% scatter(1:8,daydistance_summary(mm,:),circlesize,'r','MarkerEdgeColor','r','LineWidth',3); hold on;
% 
% end
% 
% ylim([0 500]);
% xlim([0 9]);
% 
% title('Total distance travelled rank sum test >0.05');
% xlabel('Day');
% ylabel('Distance (m)');
% 
% for rr=1:8;
%     p(rr)=ranksum(daydistance_summary(1:4,rr),daydistance_summary(5:9,rr));
% end

% group into two and plot error bars
% 
% % control group days14 or days58
% 
% control14=mean(reshape(daydistance_summary(1:4,1:4),1,16));
% control14sd=std(reshape(daydistance_summary(1:4,1:4),1,16));
% control58=mean(reshape(daydistance_summary(1:4,5:8),1,16));
% control58sd=std(reshape(daydistance_summary(1:4,5:8),1,16));
% 
% inactivation14=mean(reshape(daydistance_summary(5:9,1:4),1,20));
% inactivation14sd=std(reshape(daydistance_summary(5:9,1:4),1,20));
% inactivation58=mean(reshape(daydistance_summary(5:9,5:8),1,20));
% inactivation58sd=std(reshape(daydistance_summary(5:9,5:8),1,20));
% 
% figure;
% hold on;
% errorbar(1,control14,control14sd,'ob');
% errorbar(2,control58,control58sd,'ob');
% errorbar(1,inactivation14,inactivation14sd,'or');
% errorbar(2,inactivation58,inactivation58sd,'or');
% 
% title('Total distance travelled rank sum test >0.05');
% ylabel('Distance (m)');
% xlabel('Days');
% xlim([0.5 2.5]);
% set(gca,'XTick',[1 2]);
% set(gca,'XTickLabel',{'1-4','5-8'})


% ----------------------
% plot each box in different colours
% based on solution from http://www.mathworks.com/matlabcentral/answers/22

figure;

datacontrol=[reshape(daydistance_summary(1:4,1:4),1,16)',  reshape(daydistance_summary(1:4,5:8),1,16)',];
datacontrolposition=[1 3];
datainactivation=[reshape(daydistance_summary(5:9,1:4),1,20)', reshape(daydistance_summary(5:9,5:8),1,20)',];
datainactivationposition=[1.5 3.5];

box_control = boxplot(datacontrol,'colors','b','positions',datacontrolposition,'width',0.3); 
set(gca,'XTickLabel',{' '})

hold on;

box_inactivation = boxplot(datainactivation,'colors','r','positions',datainactivationposition,'width',0.3); 
set(gca,'XTickLabel',{' '})

xlabel('Days');
xlim([0.5 4]);
set(gca,'XTick',[1.25 3.25]);
set(gca,'XTickLabel',{'1-4','5-8'})

xlim([0.5 4]);
ylim([0 400]);

title('Mean daily distance traveled');
ylabel('Distance (m)');
xlabel('Days');


% calculate p values between each group for each day period
p14=ranksum(reshape(daydistance_summary(1:4,1:4),1,16),reshape(daydistance_summary(5:9,1:4),1,20));
p58=ranksum(reshape(daydistance_summary(1:4,5:8),1,16),reshape(daydistance_summary(5:9,5:8),1,20));


text('Position',[1.25,375],'String',sprintf('p=%s',num2str(p14,2)),'HorizontalAlignment','center');
text('Position',[3.25,375],'String',sprintf('p=%s',num2str(p58,2)),'HorizontalAlignment','center');
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/meandaydistance'),'-depsc');
close;
