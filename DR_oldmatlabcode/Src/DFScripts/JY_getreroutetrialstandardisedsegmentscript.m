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
%animals = {'M2'};
%animals = {'M2','M1','M3','K3','L2','L3'};

animals = {'N1'};

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

outall=setfilterfunction(f, 'JY_gettrajectorysegments', {'data','linpos'});
        
outall=runfilter(outall);



for mm=1:size(outall,2)
out=outall(1,mm);
uniquedays=unique(out.epochs{1,1}(:,1));
outdata=out.output{1,1};

% plot trajectory time
% setup variables
ncutoff=15; % number of traj to analyse
winsize=11;
epoch_len = [];
day_len = [];
trialsegments_all = [];
day_totalseg = [];
numseg_mean = [];
numseg_std = [];
day_totaltime = [];
totaltime_mean = cell(1,ncutoff);
totaltime_std = cell(1,ncutoff);
trials_day={};
barrier=[];
barriercount=[];




alldays=[];
% loop through days
for i = 1:size(uniquedays)
    
    jind=find(out.epochs{1,1}(:,1)==uniquedays(i,1));
    % collecting trial segments across epochs and store which trial each day and
    % epoch ends so it can be used to plot later
    tempday=[];
    for j = 1:size(jind,1);
    tempday=[tempday;out.output{1,1}{1,jind(j)}.standardisedsegments'];
    
    % collecting trial segments across epochs and store which trial each day and
    % epoch ends so it can be used to plot later
    day_len = [day_len; out.output{1,1}{1,jind(j)}.ntrial];
    barriercount=[barriercount; size(out.output{1,1}{1,jind(j)}.barrier,1)];
    barrier = [barrier; out.output{1,1}{1,jind(j)}.barrier];
   
    end
    trials_day{i}=tempday;
    alldays=[alldays;tempday];
end

% find unique trials over the entire period

    uniquetraj=unique(alldays);

    epochs = [];
    
for i=1:size(uniquedays);
    
    % go through each unique trajectoy and calculate moving average
    daytrialsummary=[];
    
    % find trajectories the 5 shortest distance
    uind=uniquetraj(1:ncutoff-1);
    for j=1:size(uind,1);
        v=uind(j);
    dayvisitones=zeros(size(trials_day{i}));
    dayvisitones((trials_day{i}==v),1)=1;
    daytrialsummary(:,j)=dayvisitones;
    end
    
    % find trajectries having 8 or more
    dayvisitones=zeros(size(trials_day{i}));
    dayvisitones((trials_day{i}>=uniquetraj(ncutoff)),1)=1;
    daytrialsummary(:,ncutoff)=dayvisitones;

    % moving average for each trajectory type [3 to 8]
    for w=1:ncutoff;
    totaltime_mean{w} = [totaltime_mean{w}; moving(daytrialsummary(:,w),winsize,@mean)];
    totaltime_std{w} = [totaltime_std{w}; moving(daytrialsummary(:,w),winsize,@std)];
    w=w+1;
    end

end





    
epochsummary=[out.epochs{1,1}(:,1) day_len];


% Plotting time per trial
figure;
set(gcf,'position',[0 0 1000 300]);
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
    text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),0.9,num2str(uniquedays(n,1)),'HorizontalAlignment','center','FontSize',12);
    set(get(get(linen,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off'); % Exclude line from legend
%     text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),1.15,[num2str(daylength(n+1,2)/100,3) 'm'],'HorizontalAlignment','center','FontSize',10);
    n=n+1;
end
hold on;
% add line to show barrier trials

barrieron=find(barrier(:,1==1));
for i=1:size(barrieron,1);
    lineb=line([barrieron(i) barrieron(i)],[0.9 0.9],'LineStyle','s', 'Color','k','LineWidth',2);
    set(get(get(lineb,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','off'); % Exclude line from legend
    hold on;
end

% plot each in a different color


col=colormap(jet(ncutoff));

%col=[255 0 0; 255 127 0;0 255 0; 0 255 255 ;255 0 255; 0 0 255];
%col=col./255;
v=1;
hold on;

for w=1:ncutoff
    
    plot(totaltime_mean{w},'Color',col(w,:),'LineWidth',3);
   
    w=w+1;
    %v=v+1;
end

%colors{8}='>=8';
%legend('3','4','5','6','7','>=8','Location','EastOutside');
%legend('Location','EastOutside');
%colorbar;
step=0.1;
colorbar('YTickLabel', num2str([uniquetraj(2):(uniquetraj(ncutoff)-uniquetraj(1))/(ceil(ncutoff/2)+1):uniquetraj(ncutoff)]',2), 'Location','EastOutside');
%colorbar;

xlim([0 max(daycumsum(:,1))]);
ylim([0 1]);
title(sprintf('Segments traversed for each trial for %s, day %d-%d \n Moving average (%d trial window)', out.animal{1,3},min(uniquedays),max(uniquedays), winsize));
xlabel('Trial');
ylabel('Proportion of moving window');

%diff=barriercount-day_len;

%clear all;
% 
% % % saving
directoryname=out.animal{1,2};
 [s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
 print(sprintf('%sPlot/behav/%s_distance',directoryname,out.animal{1,3}),'-depsc');
% 

close;
end