% calculates the number of rewards for each day


Veqn = '>=0'
minV =  str2num(Veqn(end))
maxstage = 3% [1 2 3]
minVPF = 2 %cm/sec
minPeakPF = 3
lessthan=0
includestates = 6

%Animal selection
%-----------------------------------------------------
%animals = {'K3'};

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

outall=setfilterfunction(f, 'JY_getrewards', {'data'});
        
outall=runfilter(outall);

day_velocity=[];
meanvelocity=[];
stdvelocity=[];
for mm=1:size(outall,2)
    out=outall(1,mm);
    uniquedays=unique(out.epochs{1,1}(:,1));
    outdata=out.output{1,1};

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
    
    % loop though data for each day
    for i = 1:size(uniquedays)
        trials_day = [];
        trials_day2 = [];
        jind=find(out.epochs{1,1}(:,1)==uniquedays(i,1));
        for j = 1:size(jind,1);
            trials_day=[trials_day; out.output{1,1}{1,jind(j)}.rewards];
            trials_day2=[trials_day2; out.output{1,1}{1,jind(j)}.duration];

        end

        
        daymean(mm,i)=sum(trials_day);
        daymean2(mm,i)=sum(trials_day2);
        epochs = [];
        

        
    end

end



% ----------------------
% plot each box in different colours
% based on solution from http://www.mathworks.com/matlabcentral/answers/22

daymeanvelocity=daymean./daymean2;


figure;

datacontrol=[reshape(daymeanvelocity(1:4,1:4),1,16)',  reshape(daymeanvelocity(1:4,5:8),1,16)',];
datacontrolposition=[1 3];
datainactivation=[reshape(daymeanvelocity(5:9,1:4),1,20)', reshape(daymeanvelocity(5:9,5:8),1,20)',];
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
set(gca,'YTick',[0:0.02:0.1]);
xlim([0.5 4]);
ylim([0 0.1]);

title('Rate of reward');
ylabel('Rate (rewards/second)');
xlabel('Days');


% calculate p values between each group for each day period
p14=ranksum(reshape(daymeanvelocity(1:4,1:4),1,16),reshape(daymeanvelocity(5:9,1:4),1,20));
p58=ranksum(reshape(daymeanvelocity(1:4,5:8),1,16),reshape(daymeanvelocity(5:9,5:8),1,20));


text('Position',[1.25,0.085],'String',sprintf('p=%s',num2str(p14,2)),'HorizontalAlignment','center');
text('Position',[3.25,0.085],'String',sprintf('p=%s',num2str(p58,2)),'HorizontalAlignment','center');
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/dayrewardrate'),'-depsc');
close;
