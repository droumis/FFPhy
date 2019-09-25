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
%animals = {'M2','M1','M3','K3','L2','L3','N1'};
animals = {'K3','L3','M2','N3','L2','M1','M3','N1','P2'};
%animals = {'N3'};

%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filter

days='[1:8]';%,'1:10';
%days='[1]';%,'1:10';


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

out=setfilterfunction(f, 'JY_getturnstatistic', {'data','linpos'});
        
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
    trials_day_left=[];
     
    barrier=[];
    barriercount=[];
    
    uniquedays=unique(out(1,ii).epochs{1,1}(:,1));
    
    curranimal=[];
    currraw=[];
    curranimal2={};
    % loop through days
    for i = 1:size(uniquedays)
        trials_day_left = [];
        trials_day_right=[];
        trials_day_straight=[];
        
        
        trials_duration=[];
        
        
         jind=find(out(1,ii).epochs{1,1}(:,1)==uniquedays(i,1));
%         % collecting trial segments across epochs and store which trial each day and
%         % epoch ends so it can be used to plot later
         for j = 1:size(jind,1);
%             trials_day_left=[trials_day_left; out(1,ii).output{1,1}{1,jind(j)}.leftturn];
%             trials_day_right=[trials_day_right; out(1,ii).output{1,1}{1,jind(j)}.rightturn];
%              trials_day_straight=[trials_day_straight; out(1,ii).output{1,1}{1,jind(j)}.straight];
%             
%             trials_duration=[trials_duration; out(1,ii).output{1,1}{1,jind(j)}.duration];
%             
%         end
%         % compute proportion of daily trials with [3 4 5 6 7 >8] segments
%         daytrialsummary=[];
%         daytrialraw=[];
        lengthleftturn=sum(cellfun(@(x) (x.leftturn),{out(1,ii).output{1,1}{1,jind'}}));
        lengthrightturn=sum(cellfun(@(x) (x.rightturn),{out(1,ii).output{1,1}{1,jind'}}));
        lengthstraight=sum(cellfun(@(x) (x.straight),{out(1,ii).output{1,1}{1,jind'}}));
        dayduration=sum(cellfun(@(x) (x.duration),{out(1,ii).output{1,1}{1,jind'}}));
        
%          lengthleftturn=sum(trials_day_left);
%         lengthrightturn=sum(trials_day_right);
%         lengthstraight=sum(trials_day_straight);
%         dayduration=sum(trials_duration);
        
        lengthctrleft=sum(cellfun(@(x) (x.centerturnleft),{out(1,ii).output{1,1}{1,jind'}}));
        lengthctrright=sum(cellfun(@(x) (x.centerturnright),{out(1,ii).output{1,1}{1,jind'}}));
        lengthctrstraight=sum(cellfun(@(x) (x.centerturnstraight),{out(1,ii).output{1,1}{1,jind'}}));
        lengthctrall=sum(cellfun(@(x) (x.centerturnsall),{out(1,ii).output{1,1}{1,jind'}}));
        
        lengthcnrleft=sum(cellfun(@(x) (x.cornerturnleft),{out(1,ii).output{1,1}{1,jind'}}));
        lengthcnrright=sum(cellfun(@(x) (x.cornerturnright),{out(1,ii).output{1,1}{1,jind'}}));
        lengthcnrstraight=sum(cellfun(@(x) (x.cornerturnstraight),{out(1,ii).output{1,1}{1,jind'}}));
         lengthcnrall=sum(cellfun(@(x) (x.cornerturnsall),{out(1,ii).output{1,1}{1,jind'}}));
         
    end
        
        % find trajectries having 8 or more
        daytrialsummary=lengthleftturn/(lengthleftturn+lengthrightturn);
        % daytrialraw [1.no.left 2.no.right 3.bias 4.proportion left/turns
        % 5.proportion right/turns 6.turn rate (turn/sec) 7.no.straight
        % 8.proportion left/all choices 9. proportion right/all choices 
        % 10.proportion straight/all choices;
         % 11 center left turns /center turns 
         % 12 center right turns/center turns 
         % 13 center straight/center turns 
         % 14 corner left turns /corner turns 
         % 15 corner right turns /corner turns    
         % 16 corner straight/corner turns 
         % 17 center all
         % 18 corner all
         daytrialraw=[lengthleftturn ... 
                lengthrightturn ...
            (lengthleftturn-lengthrightturn)/(lengthleftturn+lengthrightturn) ...
            lengthleftturn/(lengthleftturn+lengthrightturn) ...
            1-lengthleftturn/(lengthleftturn+lengthrightturn) ...
            (lengthleftturn+lengthrightturn)/dayduration ...
            lengthstraight  ...
            lengthleftturn/(lengthleftturn+lengthrightturn+lengthstraight) ...
        lengthrightturn/(lengthleftturn+lengthrightturn+lengthstraight) ...
         lengthstraight/(lengthleftturn+lengthrightturn+lengthstraight) ...        
         lengthctrleft/(lengthctrall) ...
         lengthctrright/(lengthctrall) ...          
         lengthctrstraight/(lengthctrall) ...        
         lengthcnrleft/(lengthcnrall) ...
         lengthcnrright/(lengthcnrall) ...
         lengthcnrstraight/(lengthcnrall)...
         lengthctrall/(lengthctrall+lengthcnrall) ...
         lengthcnrall/(lengthctrall+lengthcnrall)];
        
        
        curranimal=[curranimal;daytrialsummary];
        currraw=[currraw; daytrialraw];
        
    end
    
    animallist{ii}=curranimal;
    dayraw{ii}=currraw;
    figure;
     bar(uniquedays,dayraw{1,ii}(:,4:5),'stacked')
      title(sprintf('%s proportion left/right',out(1,ii).animal{1,3}));
%     
%     
%     %colindex=[0 51 204; 153 0 204; 0 204 153; 255 194 10; 153 204 0; 204 51 0]./255;
%     
%     colindex=[255 127 0; 152 78 163; 77 175 74; 255 255 51; 55 126 184; 228 26 28]./255;
%     
%     colormap(colindex);
%     
%     xlabel('Day')
%     ylabel('Proportion');
%     
%     
     %bar(uniquedays,dayraw{1,ii}(:,6),'stacked')
     %title(sprintf('%s trials by segment number',out(1,ii).animal{1,3}));
%     
%     legend('3','4','5','6','7','>8', 'location','EastOutside');
%     
     print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/%s_leftrightproportion',out(1,ii).animal{1,3}),'-depsc');
%     
%     close;
    
end


%% plot side by side scatter;

%% left turns


templistc=[];
for jj=1:4; % control animals
    templistc=cat(2,templistc,animallist{1,jj}(1:8,:));
end
for jj=5:9; % muscimol animals
    templistc=cat(2,templistc,animallist{1,jj}(1:8,:));
end

templistc=templistc';

datacontrol=[reshape(templistc(1:4,1:4),1,16)',  reshape(templistc(1:4,5:8),1,16)'];
datacontrolposition=[1 3];
datainactivation=[reshape(templistc(5:9,1:4),1,20)', reshape(templistc(5:9,5:8),1,20)',];
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
set(gca,'YTick',[0:0.2:1]);
xlim([0.5 4]);
ylim([0 1]);

title('Proportion of left turns');
ylabel('Proportion');
xlabel('Days');


% calculate p values between each group for each day period
p14=ranksum(reshape(daymeanvelocity(1:4,1:4),1,16),reshape(daymeanvelocity(5:9,1:4),1,20));
p58=ranksum(reshape(daymeanvelocity(1:4,5:8),1,16),reshape(daymeanvelocity(5:9,5:8),1,20));


text('Position',[1.25,18],'String',sprintf('p=%s',num2str(p14,2)),'HorizontalAlignment','center');
text('Position',[3.25,18],'String',sprintf('p=%s',num2str(p58,2)),'HorizontalAlignment','center');
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/meandayleftturnproportion'),'-depsc');
close;



%% turn bias


templistc=[];
for jj=1:4; % control animals
    templistc=cat(2,templistc,dayraw{1,jj}(1:8,3));
end
for jj=5:9; % muscimol animals
    templistc=cat(2,templistc,dayraw{1,jj}(1:8,3));
end

templistc=templistc';

datacontrol=[reshape(templistc(1:4,1:4),1,16)',  reshape(templistc(1:4,5:8),1,16)'];
datacontrolposition=[1 3];
datainactivation=[reshape(templistc(5:9,1:4),1,20)', reshape(templistc(5:9,5:8),1,20)',];
datainactivationposition=[1.5 3.5];
figure;
box_control = boxplot(datacontrol,'colors','b','positions',datacontrolposition,'width',0.3); 
set(gca,'XTickLabel',{' '})

hold on;

box_inactivation = boxplot(datainactivation,'colors','r','positions',datainactivationposition,'width',0.3); 
set(gca,'XTickLabel',{' '})

xlabel('Days');
xlim([0.5 4]);
set(gca,'XTick',[1.25 3.25]);
set(gca,'XTickLabel',{'1-4','5-8'})
set(gca,'YTick',[-1:0.5:1]);
xlim([0.5 4]);
ylim([-1 1]);

title('Bias is turn direction');
ylabel('Proportion difference');
xlabel('Days');


% calculate p values between each group for each day period
p14=ranksum(reshape(templistc(1:4,1:4),1,16),reshape(templistc(5:9,1:4),1,20));
p58=ranksum(reshape(templistc(1:4,5:8),1,16),reshape(templistc(5:9,5:8),1,20));


%text('Position',[1.25,0.9],'String',sprintf('p=%s',num2str(p14,2)),'HorizontalAlignment','center');
%text('Position',[3.25,0.9],'String',sprintf('p=%s',num2str(p58,2)),'HorizontalAlignment','center');
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnbias'),'-depsc');
close;

%% total turn rate / second

templistc=[];
for jj=1:4; % control animals
    templistc=cat(2,templistc,dayraw{1,jj}(1:8,6));
end
for jj=5:9; % muscimol animals
    templistc=cat(2,templistc,dayraw{1,jj}(1:8,6));
end

templistc=templistc';

datacontrol=[reshape(templistc(1:4,1:4),1,16)',  reshape(templistc(1:4,5:8),1,16)'];
datacontrolposition=[1 3];
datainactivation=[reshape(templistc(5:9,1:4),1,20)', reshape(templistc(5:9,5:8),1,20)',];
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
set(gca,'YTick',[0:0.1:0.5]);
xlim([0.5 4]);
ylim([0 0.2]);

title('Daily turn rate');
ylabel('Rate (turns/second)');
xlabel('Days');


% calculate p values between each group for each day period
p14=ranksum(reshape(templistc(1:4,1:4),1,16),reshape(templistc(5:9,1:4),1,20));
p58=ranksum(reshape(templistc(1:4,5:8),1,16),reshape(templistc(5:9,5:8),1,20));


text('Position',[1.25,0.18],'String',sprintf('p=%s',num2str(p14,2)),'HorizontalAlignment','center');
text('Position',[3.25,0.18],'String',sprintf('p=%s',num2str(p58,2)),'HorizontalAlignment','center');
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnrate'),'-depsc');
close;



%% all choices


% days 1-4

% get the data only for segments 4 or more, since minimal number for day
% 1-4 is 4 segments
templistc=[];
for jj=1:4; % control animals
    templistc=cat(3,templistc,dayraw{1,jj}(1:4,8:10));
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,dayraw{1,jj}(1:4,8:10));
end

% control
collectvalues=[];
collectindices=[];
for ii=1:4;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:3],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1:6],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    hold on;
end

% rearrangement indices
collectindices2=reshape([repmat([1:4]',1,3) repmat([5:8]',1,3) repmat([9:12]',1,3) repmat([13:16]',1,3)],[],1);

% lillie test for each day
resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

% do stat tests by accumarray, segments in rows

testarrayc=accumarray([collectindices collectindices2],collectvalues);
testarraycmean=mean(testarrayc,2);
testarraycstd=std(testarrayc,[],2);

 h=boxplot(collectvalues,collectindices,'colors','b','positions',[1.1:1:3.1],'width',0.2, 'symbol','o');
 set(gca,'XTickLabel',{' '})
 hold on;

%muscimol

collectvalues=[];
collectindices=[];

for ii=5:9;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:3],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1.5:1:6.5],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    hold on;
end

% rearrangement indices and rank sum test
collectindices2=reshape([repmat([1:4]',1,3) repmat([5:8]',1,3)...
    repmat([9:12]',1,3) repmat([13:16]',1,3) repmat([17:20]',1,3)],[],1);


% do stat tests by accumarray, segments in rows

testarraym=accumarray([collectindices collectindices2],collectvalues);
% lillie test for each day
resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

testarraymmean=mean(testarraym,2);
testarraymstd=std(testarraym,[],2);


h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:3.4],'width',0.2,'symbol','o');

% annotations
set(gca,'XTickLabel',{' '});

for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,1.05],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end
xlim([0.5 4]);
ylim([0 1.1]);
set(gca,'XTick',[1.25:1:3.5]);
set(gca,'XTickLabel',{'Left turn','Right turn', 'Straight'})
set(gca,'YTick',[0:0.2:1]);

title('Turn choice at intersections Day 1-4');
ylabel('Proportion');
%xlabel('Turn choice');
%figure;
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnchoice1_4'),'-depsc');
close;

combinedarraymean=[testarraycmean;testarraymmean];
combinedarraystd=[testarraycstd;testarraymstd];

%errorbar(combinedarraymean,combinedarraystd,'xb');


% days 5-8


% get the data only for segments 4 or more, since minimal number for day
% 1-4 is 4 segments
templistc=[];
for jj=1:4; % control animals
    templistc=cat(3,templistc,dayraw{1,jj}(5:8,8:10));
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,dayraw{1,jj}(5:8,8:10));
end

% control
collectvalues=[];
collectindices=[];
for ii=1:4;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:3],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1:6],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    hold on;
end

% rearrangement indices
collectindices2=reshape([repmat([1:4]',1,3) repmat([5:8]',1,3) repmat([9:12]',1,3) repmat([13:16]',1,3)],[],1);

% lillie test for each day
resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

% do stat tests by accumarray, segments in rows

testarrayc=accumarray([collectindices collectindices2],collectvalues);
testarraycmean=mean(testarrayc,2);
testarraycstd=std(testarrayc,[],2);
figure;
 h=boxplot(collectvalues,collectindices,'colors','b','positions',[1.1:1:3.1],'width',0.2, 'symbol','o');
 set(gca,'XTickLabel',{' '})
 hold on;

%muscimol

collectvalues=[];
collectindices=[];

for ii=5:9;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:3],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1.5:1:6.5],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    hold on;
end

% rearrangement indices and rank sum test
collectindices2=reshape([repmat([1:4]',1,3) repmat([5:8]',1,3)...
    repmat([9:12]',1,3) repmat([13:16]',1,3) repmat([17:20]',1,3)],[],1);


% do stat tests by accumarray, segments in rows

testarraym=accumarray([collectindices collectindices2],collectvalues);
% lillie test for each day
resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

testarraymmean=mean(testarraym,2);
testarraymstd=std(testarraym,[],2);


h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:3.4],'width',0.2,'symbol','o');


% annotations
set(gca,'XTickLabel',{' '});


for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,1.05],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end
xlim([0.5 4]);
ylim([0 1.1]);
set(gca,'XTick',[1.25:1:3.5]);

set(gca,'XTickLabel',{'Left turn','Right turn', 'Straight'});

set(gca,'YTick',[0:0.2:1]);

title('Choice at intersections Day 5-8');
ylabel('Proportion');
%xlabel('Turn choice');
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnchoice5_8'),'-depsc');
close;



combinedarraymean=[testarraycmean;testarraymmean];
combinedarraystd=[testarraycstd;testarraymstd];

%errorbar(combinedarraymean,combinedarraystd,'xb');


%% all choices by locations

catagories=6;

% days 1-4

% get the data only for segments 4 or more, since minimal number for day
% 1-4 is 4 segments
templistc=[];
for jj=1:4; % control animals
    templistc=cat(3,templistc,dayraw{1,jj}(1:4,11:16));
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,dayraw{1,jj}(1:4,11:16));
end

% control
collectvalues=[];
collectindices=[];
for ii=1:4;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1:6],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    %hold on;
end

% rearrangement indices
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories) repmat([9:12]',1,catagories) repmat([13:16]',1,catagories)],[],1);

% lillie test for each day
resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

% do stat tests by accumarray, segments in rows

testarrayc=accumarray([collectindices collectindices2],collectvalues);
testarraycmean=mean(testarrayc,2);
testarraycstd=std(testarrayc,[],2);
figure;

 h=boxplot(collectvalues,collectindices,'colors','b','positions',[1.1:1:catagories+0.1],'width',0.2, 'symbol','o');
 set(gca,'XTickLabel',{' '})
 hold on;

%muscimol

collectvalues=[];
collectindices=[];

for ii=5:9;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1.5:1:6.5],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    hold on;
end

% rearrangement indices and rank sum test
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories)...
    repmat([9:12]',1,catagories) repmat([13:16]',1,catagories) repmat([17:20]',1,catagories)],[],1);


% do stat tests by accumarray, segments in rows

testarraym=accumarray([collectindices collectindices2],collectvalues);
% lillie test for each day
resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

testarraymmean=mean(testarraym,2);
testarraymstd=std(testarraym,[],2);


h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:catagories+0.4],'width',0.2,'symbol','o');

% annotations
set(gca,'XTickLabel',{' '});

for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,1.05],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end
xlim([0.5 catagories+1]);
ylim([0 1.1]);
set(gca,'XTick',[1.25:1:catagories+0.25]);
set(gca,'XTickLabel',{'Center L','Center R', 'Center Str','Corner L','Corner R', 'Corner Str'})
set(gca,'YTick',[0:0.2:1]);

title('Turn choice at intersections Day 1-4');
ylabel('Proportion');
%xlabel('Turn choice');
%figure;
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnchoiceall1_4'),'-depsc');
close;

combinedarraymean=[testarraycmean;testarraymmean];
combinedarraystd=[testarraycstd;testarraymstd];

%errorbar(combinedarraymean,combinedarraystd,'xb');


% days 5-8


% get the data only for segments 4 or more, since minimal number for day
% 1-4 is 4 segments
templistc=[];
for jj=1:4; % control animals
    templistc=cat(3,templistc,dayraw{1,jj}(5:8,11:16));
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,dayraw{1,jj}(5:8,11:16));
end

% control
collectvalues=[];
collectindices=[];
for ii=1:4;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1:6],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    %hold on;
end

% rearrangement indices
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories) repmat([9:12]',1,catagories) repmat([13:16]',1,catagories)],[],1);

% lillie test for each day
resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

% do stat tests by accumarray, segments in rows

testarrayc=accumarray([collectindices collectindices2],collectvalues);
testarraycmean=mean(testarrayc,2);
testarraycstd=std(testarrayc,[],2);

figure;
 h=boxplot(collectvalues,collectindices,'colors','b','positions',[1.1:1:catagories+0.1],'width',0.2, 'symbol','o');
 set(gca,'XTickLabel',{' '})
 hold on;

%muscimol

collectvalues=[];
collectindices=[];

for ii=5:9;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1.5:1:6.5],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    %hold on;
end

% rearrangement indices and rank sum test
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories)...
    repmat([9:12]',1,catagories) repmat([13:16]',1,catagories) repmat([17:20]',1,catagories)],[],1);


% do stat tests by accumarray, segments in rows

testarraym=accumarray([collectindices collectindices2],collectvalues);
% lillie test for each day
resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

testarraymmean=mean(testarraym,2);
testarraymstd=std(testarraym,[],2);


h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:catagories+0.4],'width',0.2,'symbol','o');

% annotations
set(gca,'XTickLabel',{' '});

for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,1.05],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end
xlim([0.5 catagories+1]);
ylim([0 1.1]);
set(gca,'XTick',[1.25:1:catagories+0.25]);
set(gca,'XTickLabel',{'Center L','Center R', 'Center Str','Corner L','Corner R', 'Corner Str'})
set(gca,'YTick',[0:0.2:1]);

title('Turn choice at intersections Day 5-8');
ylabel('Proportion');
%xlabel('Turn choice');
%figure;
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnchoiceall5_8'),'-depsc');
close;

combinedarraymean=[testarraycmean;testarraymmean];
combinedarraystd=[testarraycstd;testarraymstd];

%errorbar(combinedarraymean,combinedarraystd,'xb');


%% turn or straight by locations

catagories=4;

% days 1-4

% get the data only for segments 4 or more, since minimal number for day
% 1-4 is 4 segments
templistc=[];
for jj=1:4; % control animals
    templistc=cat(3,templistc,[sum(dayraw{1,jj}(1:4,11:12),2) dayraw{1,jj}(1:4,13)...
        sum(dayraw{1,jj}(1:4,14:15),2) dayraw{1,jj}(1:4,16)]);
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,[sum(dayraw{1,jj}(1:4,11:12),2) dayraw{1,jj}(1:4,13)...
        sum(dayraw{1,jj}(1:4,14:15),2) dayraw{1,jj}(1:4,16)]);
end

% control
collectvalues=[];
collectindices=[];
for ii=1:4;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1:6],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    %hold on;
end

% rearrangement indices
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories) repmat([9:12]',1,catagories) repmat([13:16]',1,catagories)],[],1);

% lillie test for each day
resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

% do stat tests by accumarray, segments in rows

testarrayc=accumarray([collectindices collectindices2],collectvalues);
testarraycmean=mean(testarrayc,2);
testarraycstd=std(testarrayc,[],2);
figure;

 h=boxplot(collectvalues,collectindices,'colors','b','positions',[1.1:1:catagories+0.1],'width',0.2, 'symbol','o');
 set(gca,'XTickLabel',{' '})
 hold on;

%muscimol

collectvalues=[];
collectindices=[];

for ii=5:9;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1.5:1:6.5],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    hold on;
end

% rearrangement indices and rank sum test
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories)...
    repmat([9:12]',1,catagories) repmat([13:16]',1,catagories) repmat([17:20]',1,catagories)],[],1);


% do stat tests by accumarray, segments in rows

testarraym=accumarray([collectindices collectindices2],collectvalues);
% lillie test for each day
resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

testarraymmean=mean(testarraym,2);
testarraymstd=std(testarraym,[],2);


h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:catagories+0.4],'width',0.2,'symbol','o');

% annotations
set(gca,'XTickLabel',{' '});

for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,1.05],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end
xlim([0.5 catagories+1]);
ylim([0 1.1]);
set(gca,'XTick',[1.25:1:catagories+0.25]);
set(gca,'XTickLabel',{'Center turn','Center straight','Corner turn', 'Corner straight'})
set(gca,'YTick',[0:0.2:1]);

title('Turn choice at intersections Day 1-4');
ylabel('Proportion');
%xlabel('Turn choice');
%figure;
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnchoiceturnorstraight1_4'),'-depsc');
close;

combinedarraymean=[testarraycmean;testarraymmean];
combinedarraystd=[testarraycstd;testarraymstd];

%errorbar(combinedarraymean,combinedarraystd,'xb');


% days 5-8


% get the data only for segments 4 or more, since minimal number for day
% 1-4 is 4 segments
templistc=[];

for jj=1:4; % control animals
    templistc=cat(3,templistc,[sum(dayraw{1,jj}(5:8,11:12),2) dayraw{1,jj}(5:8,13)...
        sum(dayraw{1,jj}(5:8,14:15),2) dayraw{1,jj}(5:8,16)]);
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,[sum(dayraw{1,jj}(5:8,11:12),2) dayraw{1,jj}(5:8,13)...
        sum(dayraw{1,jj}(5:8,14:15),2) dayraw{1,jj}(5:8,16)]);
end


% control
collectvalues=[];
collectindices=[];
for ii=1:4;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1:6],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    %hold on;
end

% rearrangement indices
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories) repmat([9:12]',1,catagories) repmat([13:16]',1,catagories)],[],1);

% lillie test for each day
resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

% do stat tests by accumarray, segments in rows

testarrayc=accumarray([collectindices collectindices2],collectvalues);
testarraycmean=mean(testarrayc,2);
testarraycstd=std(testarrayc,[],2);

figure;
 h=boxplot(collectvalues,collectindices,'colors','b','positions',[1.1:1:catagories+0.1],'width',0.2, 'symbol','o');
 set(gca,'XTickLabel',{' '})
 hold on;

%muscimol

collectvalues=[];
collectindices=[];

for ii=5:9;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1.5:1:6.5],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    %hold on;
end

% rearrangement indices and rank sum test
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories)...
    repmat([9:12]',1,catagories) repmat([13:16]',1,catagories) repmat([17:20]',1,catagories)],[],1);


% do stat tests by accumarray, segments in rows

testarraym=accumarray([collectindices collectindices2],collectvalues);
% lillie test for each day
resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

testarraymmean=mean(testarraym,2);
testarraymstd=std(testarraym,[],2);


h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:catagories+0.4],'width',0.2,'symbol','o');

% annotations
set(gca,'XTickLabel',{' '});

for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,1.05],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end
xlim([0.5 catagories+1]);
ylim([0 1.1]);
set(gca,'XTick',[1.25:1:catagories+0.25]);
set(gca,'XTickLabel',{'Center turn','Center straight','Corner turn', 'Corner straight'});
set(gca,'YTick',[0:0.2:1]);

title('Turn choice at intersections Day 5-8');
ylabel('Proportion');
%xlabel('Turn choice');
%figure;
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnchoiceturnorstraight5_8'),'-depsc');
close;

combinedarraymean=[testarraycmean;testarraymmean];
combinedarraystd=[testarraycstd;testarraymstd];

%errorbar(combinedarraymean,combinedarraystd,'xb');
%% turn center or corner

catagories=2;

% days 1-4

% get the data only for segments 4 or more, since minimal number for day
% 1-4 is 4 segments
templistc=[];
for jj=1:4; % control animals
    templistc=cat(3,templistc,dayraw{1,jj}(1:4,17:18)); 
       
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,dayraw{1,jj}(1:4,17:18)); 
end

% control
collectvalues=[];
collectindices=[];
for ii=1:4;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1:6],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    %hold on;
end

% rearrangement indices
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories) repmat([9:12]',1,catagories) repmat([13:16]',1,catagories)],[],1);

% lillie test for each day
resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

% do stat tests by accumarray, segments in rows

testarrayc=accumarray([collectindices collectindices2],collectvalues);
testarraycmean=mean(testarrayc,2);
testarraycstd=std(testarrayc,[],2);
figure;

 h=boxplot(collectvalues,collectindices,'colors','b','positions',[1.1:1:catagories+0.1],'width',0.2, 'symbol','o');
 set(gca,'XTickLabel',{' '})
 hold on;

%muscimol

collectvalues=[];
collectindices=[];

for ii=5:9;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1.5:1:6.5],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    hold on;
end

% rearrangement indices and rank sum test
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories)...
    repmat([9:12]',1,catagories) repmat([13:16]',1,catagories) repmat([17:20]',1,catagories)],[],1);


% do stat tests by accumarray, segments in rows

testarraym=accumarray([collectindices collectindices2],collectvalues);
% lillie test for each day
resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

testarraymmean=mean(testarraym,2);
testarraymstd=std(testarraym,[],2);


h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:catagories+0.4],'width',0.2,'symbol','o');

% annotations
set(gca,'XTickLabel',{' '});

for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,0.75],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end
xlim([0.5 catagories+1]);
ylim([0 1]);
set(gca,'XTick',[1.25:1:catagories+0.25]);
set(gca,'XTickLabel',{'Center','Corner'})
set(gca,'YTick',[0:0.2:1]);

title('Choice by location Day 1-4');
ylabel('Proportion');
%xlabel('Turn choice');
%figure;
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnchoicelocation1_4'),'-depsc');
close;

combinedarraymean=[testarraycmean;testarraymmean];
combinedarraystd=[testarraycstd;testarraymstd];

%errorbar(combinedarraymean,combinedarraystd,'xb');


% days 5-8


% get the data only for segments 4 or more, since minimal number for day
% 1-4 is 4 segments
templistc=[];

for jj=1:4; % control animals
    templistc=cat(3,templistc,dayraw{1,jj}(5:8,17:18)); 
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,dayraw{1,jj}(5:8,17:18)); 
end


% control
collectvalues=[];
collectindices=[];
for ii=1:4;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1:6],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    %hold on;
end

% rearrangement indices
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories) repmat([9:12]',1,catagories) repmat([13:16]',1,catagories)],[],1);

% lillie test for each day
resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

% do stat tests by accumarray, segments in rows

testarrayc=accumarray([collectindices collectindices2],collectvalues);
testarraycmean=mean(testarrayc,2);
testarraycstd=std(testarrayc,[],2);

figure;
 h=boxplot(collectvalues,collectindices,'colors','b','positions',[1.1:1:catagories+0.1],'width',0.2, 'symbol','o');
 set(gca,'XTickLabel',{' '})
 hold on;

%muscimol

collectvalues=[];
collectindices=[];

for ii=5:9;
    collectvalues=[collectvalues; reshape(templistc(:,:,ii),[],1)];
    collectindices=[collectindices;reshape(flipud(repmat([1:catagories],4,1)),[],1)];
    %scatter(reshape(flipud(repmat([1.5:1:6.5],4,1)),[],1),reshape(templistc(:,:,ii),[],1));
    %hold on;
end

% rearrangement indices and rank sum test
collectindices2=reshape([repmat([1:4]',1,catagories) repmat([5:8]',1,catagories)...
    repmat([9:12]',1,catagories) repmat([13:16]',1,catagories) repmat([17:20]',1,catagories)],[],1);


% do stat tests by accumarray, segments in rows

testarraym=accumarray([collectindices collectindices2],collectvalues);
% lillie test for each day
resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

testarraymmean=mean(testarraym,2);
testarraymstd=std(testarraym,[],2);


h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:catagories+0.4],'width',0.2,'symbol','o');

% annotations
set(gca,'XTickLabel',{' '});

for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,0.85],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end
xlim([0.5 catagories+1]);
ylim([0 1]);
set(gca,'XTick',[1.25:1:catagories+0.25]);
set(gca,'XTickLabel',{'Center','Corner'});
set(gca,'YTick',[0:0.2:1]);

title('Choice by location Day 5-8');
ylabel('Proportion');
%xlabel('Turn choice');
%figure;
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnchoicelocation5_8'),'-depsc');
close;

combinedarraymean=[testarraycmean;testarraymmean];
combinedarraystd=[testarraycstd;testarraymstd];

%errorbar(combinedarraymean,combinedarraystd,'xb');

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
