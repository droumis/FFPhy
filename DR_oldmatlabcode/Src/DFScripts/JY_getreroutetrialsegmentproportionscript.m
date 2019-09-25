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
%animals = {'L2'};
%animals = {'M2','M1','M3','K3','L2','L3'};
%animals = {'M2','M1','M3','K3','L2','L3','N1','P2','N3'};

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

animallist={};
animallist2={};

for ii=1:size(out,2)
    
    % plot trajectory time
    % setup variables
    
    
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
    
    uniquedays=unique(out(1,ii).epochs{1,1}(:,1));
    
    curranimal=[];
    curranimal2={};
    % loop through days
    for i = 1:size(uniquedays)
        trials_day = [];
        jind=find(out(1,ii).epochs{1,1}(:,1)==uniquedays(i,1));
        % collecting trial segments across epochs and store which trial each day and
        % epoch ends so it can be used to plot later
        for j = 1:size(jind,1);
            trials_day=[trials_day; out(1,ii).output{1,1}{1,jind(j)}.segments];
            
        end
        % compute proportion of daily trials with [3 4 5 6 7 >8] segments
        daytrialsummary=[];
        daytriallength=size(trials_day,1);
        
        % find trajectories having 3,4,5,6,7 or more segments
        uind=[1 2 3 4 5 6 7];
        for ll=1:size(uind,2);
            v=uind(ll);
            daytrialsummary(:,v)=sum(trials_day==v)/sum(trials_day>=3);
        end
        
        
        
        % find trajectries having 8 or more
        daytrialsummary(:,8)=sum(trials_day>=8)/sum(trials_day>=3);
        
        curranimal=[curranimal;daytrialsummary];
        curranimal2{i}=trials_day;
        
    end
    
    animallist{ii}=curranimal;
    animallist2{ii}=curranimal2;
    
    
    bar(uniquedays,animallist{1,ii}(:,3:8),'stacked')
    
    
    %colindex=[0 51 204; 153 0 204; 0 204 153; 255 194 10; 153 204 0; 204 51 0]./255;
    
    colindex=[255 127 0; 152 78 163; 77 175 74; 255 255 51; 55 126 184; 228 26 28]./255;
    
    colormap(colindex);
    
    xlabel('Day')
    ylabel('Proportion');
    
    
    
    title(sprintf('%s trials by segment number',out(1,ii).animal{1,3}));
    
    legend('3','4','5','6','7','>8', 'location','EastOutside');
    
    print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/%s_segproportion',out(1,ii).animal{1,3}),'-depsc');
    
    close;
    
end


%% plot side by side scatter;

% days 1-4

% get the data only for segments 4 or more, since minimal number for day
% 1-4 is 4 segments
templistc=[];
for jj=1:4; % control animals
    templistc=cat(3,templistc,animallist{1,jj}(1:4,4:8));
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,animallist{1,jj}(1:4,4:8));
end

% control
collectvalues=[];
collectindices=[];
for ii=1:4;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:5],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1:6],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    hold on;
end

% rearrangement indices
collectindices2=reshape([repmat([1:4]',1,5) repmat([5:8]',1,5) repmat([9:12]',1,5) repmat([13:16]',1,5)],[],1);

% do stat tests by accumarray, segments in rows

testarrayc=accumarray([collectindices collectindices2],collectvalues);


% lillie test for each day
resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));


h=boxplot(collectvalues,collectindices,'colors','b','positions',[1.1:1:5.1],'width',0.2, 'symbol','o');
set(gca,'XTickLabel',{' '})
hold on;

%muscimol

collectvalues=[];
collectindices=[];

for ii=5:9;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:5],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1.5:1:6.5],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    hold on;
end

% rearrangement indices and rank sum test
collectindices2=reshape([repmat([1:4]',1,5) repmat([5:8]',1,5) repmat([9:12]',1,5) repmat([13:16]',1,5) repmat([17:20]',1,5)],[],1);


% do stat tests by accumarray, segments in rows

testarraym=accumarray([collectindices collectindices2],collectvalues);
resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

% lillie test for each day
resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:5.4],'width',0.2,'symbol','o');

% annotations
set(gca,'XTickLabel',{' '});


for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,1.05],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end
xlim([0.5 6]);
ylim([0 1.1]);
set(gca,'XTick',[1.25:1:5.5]);
set(gca,'XTickLabel',{'4','5','6','7','>8' })
set(gca,'YTick',[0 1]);



title('Proportion of daily trials with n segments Day 1-4');
ylabel('Proportion');
xlabel('Segments');


print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/proportionanalysis1_4'),'-depsc');
close;











% days 5-8
templistc=[];
for jj=1:4; % control animals
    templistc=cat(3,templistc,animallist{1,jj}(5:8,3:8));
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,animallist{1,jj}(5:8,3:8));
end

figure;

collectvalues=[];
collectindices=[];

for ii=1:4;
     collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:6],4,1)),[],1)];
    
    %scatter(reshape(flipud(repmat([1:6],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    hold on;
end

% rearrangement indices
collectindices2=reshape([repmat([1:4]',1,6) repmat([5:8]',1,6) repmat([9:12]',1,6) repmat([13:16]',1,6)],[],1);

% do stat tests by accumarray, segments in rows

testarrayc=accumarray([collectindices collectindices2],collectvalues);


% lillie test for each day
resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

h=boxplot(collectvalues,collectindices,'colors','b','positions',[1.1:1:6.1],'width',0.2,'symbol','o');
set(gca,'XTickLabel',{' '})
hold on;

collectvalues=[];
collectindices=[];

for ii=5:9;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:6],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1.5:1:6.5],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    hold on;
end

% rearrangement indices and rank sum test
collectindices2=reshape([repmat([1:4]',1,6) repmat([5:8]',1,6) repmat([9:12]',1,6) repmat([13:16]',1,6) repmat([17:20]',1,6)],[],1);


% do stat tests by accumarray, segments in rows

testarraym=accumarray([collectindices collectindices2],collectvalues);
resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:6.4],'width',0.2,'symbol','o');
set(gca,'XTickLabel',{' '})
% annotations



xlim([0.5 7]);
ylim([0 1.1]);
set(gca,'XTick',[1.25:1:6.25]);
set(gca,'XTickLabel',{'3','4','5','6','7','>8' })
set(gca,'YTick',[0 1]);

for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,1.05],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end


title('Proportion of daily trials with n segments Day 5-8');
ylabel('Proportion');
xlabel('Segments');

print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/proportionanalysis5_8'),'-depsc');
close;




% % saving
% directoryname=f.animal{1,2};
%  [s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
%  print(sprintf('%sPlot/behav/%s_bar',directoryname,out(1,ii).animal{1,3}),'-depsc');





%% old code
% save on hippo

% bar(animallist{1,1}(7,3:8))

%kstest2(animallist{1,1}(1,3:8),animallist{1,3}(1,3:8))

% h={};
% for ii=1:size(animallist2,2)
%     for jj=1:size(animallist2,2)
%         for kk=1:size(animallist2{1,jj},2)
%        h{kk}(ii,jj)=kstest2(animallist2{1,ii}{kk},animallist2{1,jj}{kk});
%         end
%     end
% end
%

% % Plotting time per trial
% figure;
% set(gcf,'position',[0 0 1000 300]);
% set(gcf,'PaperPositionMode','auto');
% set(gca,'FontSize',12);

%
%
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
% %col=colormap(hsv(6));
%
% col=[255 0 0; 255 127 0;0 255 0; 0 255 255 ;255 0 255; 0 0 255];
% col=col./255;
% v=1;
% hold on;
%
% for w=[3 4 5 6 7 8]
%
%     plot(totaltime_mean{w},'Color',col(v,:),'LineWidth',3);
%
%     w=w+1;
%     v=v+1;
% end
%
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
%
% %clear all;
%
% % % saving
% directoryname=f.animal{1,2};
%  [s,mess,messid] = mkdir(sprintf('%sPlot/behav/',directoryname));
%  print(sprintf('%sPlot/behav/%s_segments',directoryname,f.animal{1,1}),'-depsc');
%