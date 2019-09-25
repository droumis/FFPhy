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
%animals = {'P2'};
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
%epochfilter{1} = ['isequal($epoch, 6)'];
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

out=setfilterfunction(f, 'JY_intersectionanalysis', {'data','linpos'});
        
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
    
    vlengthctrleft
    vlengthctrright
    vlengthctrstraight
    
    
    
    
    
    for i = 1:size(uniquedays)
        trials_day_left = [];
        trials_day_right=[];
        trials_day_straight=[];
        
        
        trials_duration=[];
        
        
         jind=find(out(1,ii).epochs{1,1}(:,1)==uniquedays(i,1));
%         % collecting trial segments across epochs and store which trial each day and
%         % epoch ends so it can be used to plot later
         for j = 1:size(jind,1);
% 
%         lengthleftturn=sum(cellfun(@(x) (x.leftturn),{out(1,ii).output{1,1}{1,jind'}}));
%         lengthrightturn=sum(cellfun(@(x) (x.rightturn),{out(1,ii).output{1,1}{1,jind'}}));
%         lengthstraight=sum(cellfun(@(x) (x.straight),{out(1,ii).output{1,1}{1,jind'}}));
%         dayduration=sum(cellfun(@(x) (x.duration),{out(1,ii).output{1,1}{1,jind'}}));
%         

        
        lengthctrleft=mean(cellfun(@(x) (x.centerleftvelocity),{out(1,ii).output{1,1}{1,jind'}}));
        lengthctrright=mean(cellfun(@(x) (x.centerrightvelocity),{out(1,ii).output{1,1}{1,jind'}}));
        lengthctrstraight=mean(cellfun(@(x) (x.centerstraightvelocity),{out(1,ii).output{1,1}{1,jind'}}));
        %lengthctrall=sum(cellfun(@(x) (x.centersall),{out(1,ii).output{1,1}{1,jind'}}));
        
        lengthcnrleft=mean(cellfun(@(x) (x.cornerleftvelocity),{out(1,ii).output{1,1}{1,jind'}}));
        lengthcnrright=mean(cellfun(@(x) (x.cornerrightvelocity),{out(1,ii).output{1,1}{1,jind'}}));
        lengthcnrstraight=mean(cellfun(@(x) (x.cornerstraightvelocity),{out(1,ii).output{1,1}{1,jind'}}));
         %lengthcnrall=sum(cellfun(@(x) (x.cornersall),{out(1,ii).output{1,1}{1,jind'}}));
% 
        vlengthctrleft=[vlengthctrleft; cellfun(@(x) (x.centerleftvelocity),{out(1,ii).output{1,1}{1,jind'}})];
        vlengthctrright=[vlengthctrright; cellfun(@(x) (x.centerrightvelocity),{out(1,ii).output{1,1}{1,jind'}})];
        vlengthctrstraight=[vlengthctrstraight; cellfun(@(x) (x.centerstraightvelocity),{out(1,ii).output{1,1}{1,jind'}})];
        
        vlengthcnrleft=[vlengthcnrleft;cellfun(@(x) (x.cornerleftvelocity),{out(1,ii).output{1,1}{1,jind'}})];
        vlengthcnrright=[vlengthcnrright;cellfun(@(x) (x.cornerrightvelocity),{out(1,ii).output{1,1}{1,jind'}})];
        vlengthcnrstraight=[vlengthcnrstraight;cellfun(@(x) (x.cornerstraightvelocity),{out(1,ii).output{1,1}{1,jind'}})];


         poscornerindex=cell2mat(cellfun(@(x) (x.choicestartendindx(x.cornerindex,:)),...
             {out(1,ii).output{1,1}{1,jind'}},'UniformOutput', false)');
         poscenterindex=cell2mat(cellfun(@(x) (x.choicestartendindx(x.centerindex,:)),...
             {out(1,ii).output{1,1}{1,jind'}},'UniformOutput', false)');
         
         
         pos=cell2mat(cellfun(@(x) x.pos(:,2:3),{out(1,ii).output{1,1}{1,jind'}},'UniformOutput', false)');
         
%          for iii=1:50;
%              
%              plot(pos(poscornerindex(iii,1):poscornerindex(iii,2),1),...
%                  pos(poscornerindex(iii,1):poscornerindex(iii,2),2),'b')
%              hold on;
%          end
         
         for iii=1:7;
             
             plot(pos(poscenterindex(iii,1):poscenterindex(iii,2),1),...
                 pos(poscenterindex(iii,1):poscenterindex(iii,2),2),'r')
             hold on;
         end
    end
        
        % find trajectries having 8 or more
        %daytrialsummary=lengthleftturn/(lengthleftturn+lengthrightturn);
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
%          daytrialraw=[lengthleftturn ... 
%                 lengthrightturn ...
%             (lengthleftturn-lengthrightturn)/(lengthleftturn+lengthrightturn) ...
%             lengthleftturn/(lengthleftturn+lengthrightturn) ...
%             1-lengthleftturn/(lengthleftturn+lengthrightturn) ...
%             (lengthleftturn+lengthrightturn)/dayduration ...
%             lengthstraight  ...
%             lengthleftturn/(lengthleftturn+lengthrightturn+lengthstraight) ...
%         lengthrightturn/(lengthleftturn+lengthrightturn+lengthstraight) ...
%          lengthstraight/(lengthleftturn+lengthrightturn+lengthstraight) ...        
%          lengthctrleft/(lengthctrall) ...
%          lengthctrright/(lengthctrall) ...          
%          lengthctrstraight/(lengthctrall) ...        
%          lengthcnrleft/(lengthcnrall) ...
%          lengthcnrright/(lengthcnrall) ...
%          lengthcnrstraight/(lengthcnrall)...
%          lengthctrall/(lengthctrall+lengthcnrall) ...
%          lengthcnrall/(lengthctrall+lengthcnrall)];
           
     
     daytrialraw=[...      
         lengthctrleft ...
         lengthctrright ...          
         lengthctrstraight ...        
         lengthcnrleft...
         lengthcnrright ...
         lengthcnrstraight];
  
        
        %curranimal=[curranimal;daytrialsummary];
        currraw=[currraw; daytrialraw];
        
    end
    
    %animallist{ii}=curranimal;
    dayraw{ii}=currraw;
%     figure;
%      bar(uniquedays,dayraw{1,ii}(:,4:5),'stacked')
%       title(sprintf('%s trials by segment number',out(1,ii).animal{1,3}));
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
%     %print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/%s_segproportion',out(1,ii).animal{1,3}),'-depsc');
%     
%     close;
    
end



%% Plotting all choices




%% all choices by locations

catagories=6;

% days 1-4

% get the data only for segments 4 or more, since minimal number for day
% 1-4 is 4 segments
templistc=[];
for jj=1:4; % control animals
    templistc=cat(3,templistc,dayraw{1,jj}(1:4,1:6));
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,dayraw{1,jj}(1:4,1:6));
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
%resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

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
%resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

testarraymmean=mean(testarraym,2);
testarraymstd=std(testarraym,[],2);


h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:catagories+0.4],'width',0.2,'symbol','o');

% annotations
set(gca,'XTickLabel',{' '});

for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,68.05],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end
xlim([0.5 catagories+1]);
ylim([0 70]);
set(gca,'XTick',[1.25:1:catagories+0.25]);
set(gca,'XTickLabel',{'Center L','Center R', 'Center Str','Corner L','Corner R', 'Corner Str'})
set(gca,'YTick',[0:10:70]);

title('Mean turn speed at intersections Day 1-4');
ylabel('Speed cm/s');
%xlabel('Turn choice');
%figure;
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnspeed1_4'),'-depsc');
%close;

combinedarraymean=[testarraycmean;testarraymmean];
combinedarraystd=[testarraycstd;testarraymstd];

%errorbar(combinedarraymean,combinedarraystd,'xb');


% days 5-8


% get the data only for segments 4 or more, since minimal number for day
% 1-4 is 4 segments
templistc=[];
for jj=1:4; % control animals
    templistc=cat(3,templistc,dayraw{1,jj}(5:8,1:6));
end
for jj=5:9; % muscimol animals
    templistc=cat(3,templistc,dayraw{1,jj}(5:8,1:6));
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
%resultsc=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

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
%resultsm=accumarray(collectindices,collectvalues,[],@(x) lillietest(x));

testarraymmean=mean(testarraym,2);
testarraymstd=std(testarraym,[],2);


h=boxplot(collectvalues,collectindices,'colors','r','positions',[1.4:1:catagories+0.4],'width',0.2,'symbol','o');

% annotations
set(gca,'XTickLabel',{' '});

for ii=1:size(testarrayc,1);
    ranksumresults(ii)=ranksum(testarrayc(ii,:),testarraym(ii,:));

text('Position',[ii+0.25,68.05],'String',sprintf('p=%s',num2str(ranksumresults(ii),2)),'HorizontalAlignment','center');

    

end
xlim([0.5 catagories+1]);
ylim([0 70]);
set(gca,'XTick',[1.25:1:catagories+0.25]);
set(gca,'XTickLabel',{'Center L','Center R', 'Center Str','Corner L','Corner R', 'Corner Str'})
set(gca,'YTick',[0:10:70]);

title('Mean turn speed at intersections Day 5-8');
ylabel('Speed cm/s');
%xlabel('Turn choice');
%figure;
print(sprintf('/home/jai/Documents/Projects/DecisionMaking/ACCInactivation/turnspeed5_8'),'-depsc');
%close;

combinedarraymean=[testarraycmean;testarraymmean];
combinedarraystd=[testarraycstd;testarraymstd];

%errorbar(combinedarraymean,combinedarraystd,'xb');

