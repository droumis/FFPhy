% calculates segment transition at the first intersection of each trial


Veqn = '>=0'
minV =  str2num(Veqn(end))
maxstage = 3% [1 2 3]
minVPF = 2 %cm/sec
minPeakPF = 3
lessthan=0
includestates = 6

%Animal selection
%-----------------------------------------------------
animals = {'CML21'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filter

days='[09]';%,'1:10';



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

out=setfilterfunction(f, 'JY_gettrajectorysegments', {'data','linpos'});
        
out=runfilter(out);


% plot trajectory time
% setup variables
winsize=3;
epoch_len = [];
day_len = [];
trialsegments_all = [];
day_totalseg = [];
numseg_mean = [];
numseg_std = [];
day_totaltime = [];
totaltime_mean = cell(1,8);
totaltime_std = cell(1,8);
trials_day=[];
barrier=[];
barriercount=[];

uniquedays=unique(out.epochs{1,1}(:,1));

% loop through days
dayvisitones=[];
for i = 1:size(uniquedays)
    trials_day = [];
    jind=find(out.epochs{1,1}(:,1)==uniquedays(i,1));
    % collecting trial segments across epochs and store which trial each day and
    % epoch ends so it can be used to plot later
    for j = 1:size(jind,1);
    trials_day=[trials_day; out.output{1,1}{1,jind(j)}.segment12];
    % collecting trial segments across epochs and store which trial each day and
    % epoch ends so it can be used to plot later
    day_len = [day_len; out.output{1,1}{1,jind(j)}.ntrial];
    barriercount=[barriercount; size(out.output{1,1}{1,jind(j)}.barrier,1)];
    barrier = [barrier; out.output{1,1}{1,jind(j)}.barrier];
   
    end

    epochs = [];
    
    % find unique 1st segment trajectories
    
    unique1st=unique(trials_day(:,1));
    unique2ndout={};
    % go through each unique trajectoy and calculate moving average
    daytrialsummary={};
    segment1dir1=[];
    segment1dir2=[];
    
    % find the 2 starting segments
    for i=1:size(unique1st,1);
        % find the 2nd segment 
        daytrialsummary{i}=[];
        unique2nd=unique(trials_day(trials_day(:,1)==unique1st(i),2));
        unique2ndout{i}=unique2nd;
        for j=1:size(unique2nd,1)
            dayvisitones=zeros(size(trials_day(trials_day(:,1)==unique1st(i),2),1),1);
            dayvisitones(trials_day(trials_day(:,1)==unique1st(i),2)==unique2nd(j))=1;
            
            daytrialsummary{i}(:,j)=dayvisitones;
        
        end
    end
    
    
    for w=1:size(unique1st,1)
    totaltime_mean{w} = [totaltime_mean{w}; moving(daytrialsummary{w},winsize,@mean)];
    totaltime_std{w} = [totaltime_std{w}; moving(daytrialsummary{w},winsize,@std)];
    w=w+1;
    end
end
    
epochsummary=[out.epochs{1,1}(:,1) day_len];


% Plotting time per trial
figure;
set(gcf,'position',[0 0 800 300]);
set(gcf,'PaperPositionMode','auto');
set(gca,'FontSize',12);

% make grey/white strips
% o=1;
% col=[];

% uniquedays=unique(epochsummary(:,1));
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
%     text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),0.9,num2str(uniquedays(n,1)),'HorizontalAlignment','center','FontSize',12);
%     set(get(get(linen,'Annotation'),'LegendInformation'),...
%         'IconDisplayStyle','off'); % Exclude line from legend
% %     text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),1.15,[num2str(daylength(n+1,2)/100,3) 'm'],'HorizontalAlignment','center','FontSize',10);
%     n=n+1;
% end
% hold on;
% % add line to show barrier trials
% 
% barrieron=find(barrier(:,1==1));
% for i=1:size(barrieron,1);
%     lineb=line([barrieron(i) barrieron(i)],[0.9 0.9],'LineStyle','s', 'Color','k','LineWidth',2);
%     set(get(get(lineb,'Annotation'),'LegendInformation'),...
%         'IconDisplayStyle','off'); % Exclude line from legend
%     hold on;
% end
% 
% % plot each in a different color
% 
col=colormap(lines(10));

% col=[208 52 39;6 95 133;33 148 89;50 82 0;181 79 36; 79 16 199];
% col=col./255;
v=1;
hold on;

for w=[1]
    for x=1:size(totaltime_mean{w},2)
    
    plot(totaltime_mean{w}(:,x),'Color',col(v,:),'LineWidth',3);
    v=v+1;
    end  
    w=w+1;
    
    
end
tmpStr=sprintf('segment %g;', unique2ndout{1});
C = dataread('string', tmpStr, '%s', 'delimiter', ';'); 

legend(C,'Location','EastOutside');




title(sprintf('1st segment choice departing from segment %s for %s day %s \n Moving average (%s trial window)',...
    num2str(unique1st(1)),out.animal{1,1}, num2str(days), num2str(winsize)));





figure;
set(gcf,'position',[0 0 800 300]);
set(gcf,'PaperPositionMode','auto');
set(gca,'FontSize',12);

hold on;

for w=[2]
    for x=1:size(totaltime_mean{w},2);
        
    plot(totaltime_mean{w}(:,x),'--','Color',col(x,:),'LineWidth',3);
    v=v+1;
    end  
    w=w+1;
    
    
end

title(sprintf('1st segment choice departing from segment %s for %s day %s \n Moving average (%s trial window)',...
    num2str(unique1st(2)),out.animal{1,1}, num2str(days), num2str(winsize)));

tmpStr=sprintf('segment %g;', unique2ndout{2});
C = dataread('string', tmpStr, '%s', 'delimiter', ';'); 

legend(C,'Location','EastOutside');

% colors{8}='>=8';
% legend('3','4','5','6','7','>=8','Location','EastOutside');




% 
% xlim([0 max(daycumsum(:,1))]);
% ylim([0 1]);
% title(sprintf('Segments traversed for each trial for %s, day %d-%d \n Moving average (%d trial window)', out.animal{1,1},min(uniquedays),max(uniquedays), winsize));
% xlabel('Trial');
% ylabel('Proportion of moving window');
% 
% diff=barriercount-day_len;

clear all;

% % saving
% [s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
% print(sprintf('%sPlot/behav/%s_trialtimeavg_d%d-%d_w%d',directoryname,fileprefix,days(1),days(end),winsize),'-depsc');

