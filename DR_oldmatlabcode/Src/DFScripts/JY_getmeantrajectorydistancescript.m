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

%animals = {'K3','L3','M2','N1','L2','M3','M1'};
animals = {'K3','L3','M2','N3','L2','M1','M3','N1','P2'};

%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filter

days='[1:8]';%,'1:10';

%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = ['isequal($epochtype, ''Run'')'];
epochfilter{1} = ['ismember($epoch, [2 4 6])'];
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
% 
% uniquedays=unique(out.epochs{1,1}(:,1));
% daydistance=[];
% loop though data for each day

% get the mean distance
meandistance=[];
meantrialsegments=[];
rawcontrol=[];
rawinactivation=[];


%% mean distance for each epoch

trajectorydisplacement=cell(1,24);

for ii=1:4
    
    for jj=1:24
        
        trajectorydisplacement{1,jj}=[trajectorydisplacement{1,jj};...
            out(1,ii).output{1,1}{1,jj}.trajectorydisplacement];
    end
    %meandistance(ii,:)=cell2mat(cellfun(@(x) x.meanepochtrajectorydisplacement, out(1,ii).output{1,1},'UniformOutput', false));
    
    %meantrialsegments(ii,:)=cell2mat(cellfun(@(x) x.meantrialsegments, out(1,ii).output{1,1},'UniformOutput', false));
end

controlm=cellfun(@(x) mean(x),trajectorydisplacement);
controls=cellfun(@(x) std(x)/sqrt(size(x,1)),trajectorydisplacement);


trajectorydisplacement=cell(1,24);
for ii=5:9
    
    for jj=1:24
        
        trajectorydisplacement{1,jj}=[trajectorydisplacement{1,jj};...
            out(1,ii).output{1,1}{1,jj}.trajectorydisplacement];
    end
    %meandistance(ii,:)=cell2mat(cellfun(@(x) x.meanepochtrajectorydisplacement, out(1,ii).output{1,1},'UniformOutput', false));
    
    %meantrialsegments(ii,:)=cell2mat(cellfun(@(x) x.meantrialsegments, out(1,ii).output{1,1},'UniformOutput', false));
end
inactivationm=cellfun(@(x) mean(x),trajectorydisplacement);
inactivations=cellfun(@(x) std(x)/sqrt(size(x,1)),trajectorydisplacement);

figure;
errorbar([1:24],controlm,controls,'ob');
hold on;
errorbar([1:24],inactivationm,inactivations,'or');


xlim([0 25]);
ylim([0 1400]);
set(gca,'XTick',[2:3:25]);
set(gca,'XTickLabel',{'1','2', '3','4','5', '6', '7', '8'})


title('Mean epoch trajectory distance Day 1-8 +/-sem');
ylabel('Distance cm');
xlabel('Day');
%figure;
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/meanepochdistance'),'-depsc');
%close;





%control
% 
% control=mean(meandistance(1:4,:),1);
% inactivation=mean(meandistance(5:9,:),1);
% 
% plot(control,'b');
% scatter(reshape(repmat(1:24,4,1),1,[]),reshape(meandistance(1:4,:),1,[]),'b');
% hold on;
% scatter(reshape(repmat(1:24,5,1),1,[]),reshape(meandistance(5:9,:),1,[]),'r');
% ylim([0 700]);
% 
% hold on;
% plot(inactivation,'r');
% bar(inactivation,'r');

 
% Plotting time per trial
% figure;
% set(gcf,'position',[0 0 800 400]);
% set(gcf,'PaperPositionMode','auto');
% set(gca,'FontSize',12);
% 
% figure;
% plot(meandistance(1:4,:)','or');
% 
% hold on;
% plot(meandistance(5:7,:)','ob');
% % 
% % 
% % figure;
% % plot(mean(meandistance(1:3,:),1)','-r');
% 
% figure;
% errorbar(nanmean(meandistance(1:4,:),1)',nanstd(meandistance(1:4,:),1)','or');
% hold on;
% errorbar(nanmean(meandistance(5:7,:),1)',nanstd(meandistance(5:7,:),1)','ob');
% 
% figure;
% errorbar(nanmean(meantrialsegments(1:4,:),1)',nanstd(meantrialsegments(1:4,:),1)','or');
% hold on;
% errorbar(nanmean(meantrialsegments(5:7,:),1)',nanstd(meantrialsegments(5:7,:),1)','ob');

% annotating

% for ii = 1:length(day_len)
%     plot([sum(day_len(1:ii)) sum(day_len(1:ii))]-(winsize-1)*(ii), [0 200], 'k');
%     %text(sum(day_len(1:ii))-(winsize)*ii, 85, sprintf('Day %d\nn=%d\ntime=%d',days(ii),day_len(ii),floor(day_totaltime(ii))), ...
%         %'HorizontalAlignment','right','FontSize',10);
%     
% end


% title(sprintf('Trial distance for %s, day %d-%d \n Moving average (%d trial window)', out.animal{1,1},min(uniquedays),max(uniquedays), winsize));
% xlabel('Trial');
% ylabel('Multiples of shortest distance');
% 
% [s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
% print(sprintf('%sPlot/behav/%s_distance',directoryname,f.animal{1,1}),'-depsc');
% 
% 
figure;

bar(daydistance./100);

title(sprintf('Distance travelled for %s, day %d-%d', out.animal{1,1},min(uniquedays),max(uniquedays)));
xlabel('Day');
ylabel('Distance (m)');
% 
% %clear all;
% 
% % % saving
% 
% directoryname=f.animal{1,2};
% [s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
% print(sprintf('%sPlot/behav/%s_totaldistance',directoryname,f.animal{1,1}),'-depsc');