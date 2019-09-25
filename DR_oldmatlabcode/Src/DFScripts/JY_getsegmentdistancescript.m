% shows the linear displacement (trajectory distance) of each trial grouped by number of segments
% taken

% shows the 2D displacement (trajectory displacement) of each trial grouped by number of segments
% taken


Veqn = '>=0'
minV =  str2num(Veqn(end))
maxstage = 3% [1 2 3]
minVPF = 2 %cm/sec
minPeakPF = 3
lessthan=0
includestates = 6

%Animal selection
%-----------------------------------------------------
animals = {'M1'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filter

days='[1:10]';%,'1:10';

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

out=setfilterfunction(f, 'JY_gettrajectorydistance', {'data','linpos'});

out=runfilter(out);


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



% loop though data for each day
trials_day = [];
trials_day2 = [];
trials_day3= {};
trials_day_ad = [];
trials_day2_ad = [];
trials_day3_ad= {};


for i = 1:size(out.epochs{1,1}(:,1),1)
    
    
    epochindex=min(size(out.output{1,1}{1,i}.trialsegments,1), size(out.output{1,1}{1,i}.trajectorydistance,1));
    epochindexad=min(size(out.output{1,1}{1,i}.trialsegments,1), size(out.output{1,1}{1,i}.trajectorydisplacement,1));
    trials_day=[trials_day; out.output{1,1}{1,i}.trialsegments(1:epochindex,1)];
    trials_day_ad=[trials_day_ad; out.output{1,1}{1,i}.trialsegments(1:epochindexad,1)];
    
    
    trials_day2=[trials_day2; out.output{1,1}{1,i}.trajectorydistance(1:epochindex,1)];
    trials_day2_ad=[trials_day2_ad; out.output{1,1}{1,i}.trajectorydisplacement(1:epochindexad,1)];
    trials_day3=[trials_day3 out.output{1,1}{1,i}.trialsegments_details(1, 1:epochindex)];
    trials_day3_ad=[trials_day3_ad out.output{1,1}{1,i}.trialsegments_details(1, 1:epochindexad)];
    
    % collecting trial segments across epochs and store which trial each day and
    % epoch ends so it can be used to plot later
    
    
    
    % Collect the time of each trial across all epochs in the day
    
    %     totaltime_mean = [totaltime_mean movingstat(trial_time',winsize,@mean)'];
    %     totaltime_std = [totaltime_std movingstat(trial_time',winsize,@std)'];
    
    %trialsegments_all = [trialsegments_all trialsegments_day];
    
end


% find unique segment numbers
uniqueindex=unique(trials_day);
dsegs={};
for i=1:size(uniqueindex,1)
    dtable{i}=trials_day2(trials_day==uniqueindex(i));
    dtablead{i}=trials_day2_ad(trials_day_ad==uniqueindex(i));
end



dstat(1,:)=cellfun(@(x) mean(x,1),dtable);
dstat(2,:)=cellfun(@(x) std(x,1)/sqrt(size(x,1)),dtable);
dstat(3,:)=cellfun(@(x) size(x,1),dtable);

dstatad(1,:)=cellfun(@(x) mean(x,1),dtablead);
dstatad(2,:)=cellfun(@(x) std(x,1)/sqrt(size(x,1)),dtablead);
dstatad(3,:)=cellfun(@(x) size(x,1),dtablead);

% Plotting time per trial
figure;
set(gcf,'position',[0 0 800 400]);
set(gcf,'PaperPositionMode','auto');
set(gca,'FontSize',12);

% make grey/white strips

%hist(dtable{1,1},100)
%bar(dstat(3,:))
boxplot(trials_day2, trials_day)
xlabel('Trials with n segments');
ylabel('cm');
title(sprintf('Distribution of trial linear distances for grouped by trial segments \n %s day %s to %s',...
    animals{1,1},num2str(min(str2num(days))),num2str(max(str2num(days)))));
%,uniqueindex)
[p_KW,table,stats]=kruskalwallis(trials_day2, trials_day);
[p_RS,table,stats]=ranksum(dtable{1,1}(:),dtable{1,2}(:));

figure;
% 2D displacement
boxplot(trials_day2_ad, trials_day_ad)
xlabel('Trials with n segments');
ylabel('cm');
title(sprintf('Distribution of trial 2D displacement for grouped by trial segments \n %s day %s to %s',...
    animals{1,1},num2str(min(str2num(days))),num2str(max(str2num(days)))));
%,uniqueindex)
[adp_KW,table,adstats]=kruskalwallis(trials_day2_ad, trials_day_ad);
[adp_RS,adtable,adstats]=ranksum(dtablead{1,1}(:),dtablead{1,2}(:));




figure;
bar(dstat(3,:))
%set(gca,'XTick',[])
set(gca,'XTickLabel',mat2cell(uniqueindex',1,length(uniqueindex)))
title(sprintf('Distribution of trial distances for grouped by trial segments \n %s day %s to %s',...
    animals{1,1},num2str(min(str2num(days))),num2str(max(str2num(days)))));
xlabel('Trials n segments');
ylabel('Counter');

clear all;

%errorbar(dstat(1,:),dstat(2,:),'x')

% find the first segment choice for each 



