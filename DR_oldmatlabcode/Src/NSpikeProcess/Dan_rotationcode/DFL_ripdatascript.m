% uses 2 algorithms to plot reactivation probability

% anabelle's getripactivprob
% JY's own way of calculating


global subplot_count;
subplot_count = 1;

Veqn = '>=0';
minV =  str2num(Veqn(end));
maxstage = 3; % [1 2 3]
minVPF = 3; %cm/sec
minPeakPF = 3;
lessthan=0;
includestates = 6;

%Animal selection
%-----------------------------------------------------
animals = {'I1'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filterionno

days='[15]';%,'1:10';
%days = '[1:1]';
%days = '[9:9]';

%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

%epochfilter{1} = [''];

epochfilter = ['isequal($epochtype, ''Run'')'];
%epochfilter = ['isequal($epoch, 6)'];

%cellfilter = '((isequal($area, ''CA1'') && ($spikewidth <10))'; % && ($data>50)'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
cellfilter = '(isequal($area, ''CA1'') && ($meanrate <20 ) && ($numspikes>50) && ($spikewidth<15))'  ;

%timefilter = { {'JY_getlinvelocity', '($velocity <0.07) & ($velocity >-0.07)'} };
timefilter = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
   {'JY_getlinvelocity', strcat('abs($velocity < ',num2str(minVPF),')')}...
   {'JY_getbarrier','($barrier== 0)'}};

% timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%     {'JY_getlinvelocity', strcat('abs($velocity > ',num2str(minVPF),')')}...
%     {'JY_getbarrier','($barrier== 0)'}};

timefilter2 = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
   {'JY_getlinvelocity', strcat('abs($velocity < ',num2str(minVPF),')')}...
   {'JY_getbarrier','($barrier== 1)'}};

% timefilter2 = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%     {'JY_getlinvelocity', strcat('abs($velocity > ',num2str(minVPF),')')}...
%     {'JY_getbarrier','($barrier== 1)'}};

%timefilter3 = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},...
%    {'JY_getlinvelocity', strcat('abs($velocity < ',num2str(minVPF),')')}};


%timefilter = { {'JY_getriptimes','($nripples > 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { {'JY_getlinvelocity', '(($velocity) >= 0))', 6} };
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
f2 = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter2);
%f3 = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter3);


%only include cells with placefields
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator3 = 'multicellanal';
%iterator3 = 'JY_singlecellanal';

%iterator = 'JY_singleepochcellanal';
%iterator = 'multicellanal';
%iterator = 'JY_singleepochanal';
iterator = 'JY_singlecellanal';

f = setfilteriterator(f,iterator);
f2 = setfilteriterator(f2,iterator);
%f3 = setfilteriterator(f3,iterator3);


% out = plottrajdata(index, excludetimes, spikes, linpos, includestates, minV, varargin)

f = setfilterfunction(f, 'getripactivprob', {'spikes','ripples','task','cellinfo','cellfilter'});
f2 = setfilterfunction(f2, 'getripactivprob', {'spikes','ripples','task','cellinfo','cellfilter'});
%f3 = setfilterfunction(f3, 'DFL_ripdata', {'spikes','data'}, includestates, minV, f3);
%f4=setfilterfunction(f, 'getriprate', {'ripples'});

f = runfilter(f);
f2 = runfilter(f2);
%f4 = runfilter(f4);

%f3 = runfilter(f3);


% barrier cells for I1

barriercelllist=[10 8 3;10 11 3; 10 12 1; 10 13 6; 12 5 4; 12 5 7; 12 6 5; 12 7 3; 12 8 2; 12 8 4;...
                 13 6 3; 13 6 6; 13 8 2; 13 8 4; 13 13 3; 14 5 5; 15 5 1];
             


%% Compare barrier vs no barrier
%plot and loop through indices
list=[2:2:2*size(f.epochs{1,1},1)];
boxgroup=[];

% define which cells are barrier cells

barriercells=barriercelllist(find(barriercelllist(:,1)==unique(f.epochs{1,1}(:,1))),2:3);

% give these an index of 1

for i=1:size(f.epochs{1,1},1)
    index=zeros(size(f.data{1,1}{1,i},1),1);
    index(rowfind(f.data{1,1}{1,i},barriercells)~=0)=1;
    index=index+1;
    
    boxgroup=[boxgroup;index];
end





for i=1:size(f2.epochs{1,1},1)
    
    index=zeros(size(f2.data{1,1}{1,i},1),1);
    index(rowfind(f2.data{1,1}{1,i},barriercells)~=0)=1;
    index=index+3;
    
    
    boxgroup=[boxgroup;index];
end

% 
 day=f.epochs{1,1}(1,1);
% epoch=f.epochs{1,1}(1,2);

boxdata=[f.output{1,1}(:,1);f2.output{1,1}(:,1)];
%boxgroup=[zeros(size(f.output{1,1}(:,1)));ones(size(f2.output{1,1}(:,1)))];

figure;
%boxplot(boxdata,boxgroup,'labels',{'no barrier','barrier'});



%boxplot(boxdata,boxgroup,'labels',{'Ep 2 -','Ep 2 +','Ep 4 -','Ep 4 +','Ep 6 -','Ep 6 +'});

boxplot(boxdata,boxgroup,'labels',{'NB-NB','NB-B','B-NB','B-B'});


%xlabel('+/- Barrier');
ylabel('Reactivation probability');


title(sprintf('Ripple activation probability for %s \n day %s ',animals{1,1},num2str(day)));

% ----Saving----
    % Saves figure as pdf
    % First checks if a folder called Plot exists in the processed data folder,
    % if not, creates it.
    
    cd(f.animal{1,2});
    plotdir = dir('Plot');
    if (isempty(plotdir))
        %an a plot folder needs to be created
        !mkdir Plot
    end
    
    % change to that directory and saves the figure with file name
    % animal_day_epoch
    cd(strcat(f.animal{1,2},'Plot/'));
    figurename = strcat(animals{1,1},'_ripactivprob_NB_B',num2str(day));
    
    saveas(gcf, figurename, 'fig');
    
    %Closes the figure
    close;

    
%     %% no barrier data
%     
%     % ger day indices
%     days=unique(f.epochs{1,1}(:,1));
%     % get number of cells per epoch
%     datanumber=cellfun(@(x) size(x,1),f.data{1,1});
%     genindex=[];
%     % generate an index for values
%     for i=1:size(f.epochs{1,1})
%         tmpind=ones(datanumber(i),1)*f.epochs{1,1}(i,1);
%         genindex=[genindex;tmpind];
%     end
%         
%     
%     
% % get average ripple reactivation rate for each day
% 
% boxplot(f.output{1,1},genindex);
% xlabel('Day');
% ylabel('Reactivation probability');
% title(sprintf('Ripple activation probability for %s \n day %s at velocity < %s cm/s ',animals{1,1},num2str(day),num2str(minVPF)));
%     
%     
% %% riprate
% 
% % ger day indices
%     days=unique(f4.epochs{1,1}(:,1));
%     % get number of cells per epoch
%     datanumber=cellfun(@(x) size(x,1),f4.data{1,1});
%     genindex=[];
%     % generate an index for values
%     for i=1:size(f4.epochs{1,1})
%         tmpind=ones(datanumber(i),1)*f4.epochs{1,1}(i,1);
%         genindex=[genindex;tmpind];
%     end
%         
%     
%     
% % get average ripple reactivation rate for each day
% 
% boxplot(f4.output{1,1}(:,1),genindex);
% xlabel('Day');
% ylabel('Ripple rate');
% title(sprintf('Ripple rate for %s \n',animals{1,1}));
%     
% 
% 
%     
% clear all;

% 
% if isequal(iterator3,'multicellanal');
%  results=[];   
% for i=1:size(f3.epochs{1,1},1)
%     % find all non-barrier
%     
%     nbspikesind=find(f3.output{1}(i).ripple_barrierlist==0);
%     nbspikes=f3.output{1}(i).ripple_cell_count(nbspikesind);
%     nbreactivation=sum(nbspikes>0)/size(nbspikes,2);
%     % find all barrier
%     bspikesind=find(f3.output{1}(i).ripple_barrierlist==1);
%     breactivation=0;
%     bspikes=[];
%     if ~isempty(bspikesind)
%         bspikes=f.output{1}(i).ripple_cell_count(bspikesind);
%         breactivation=sum(bspikes>0)/size(bspikes,2);
%     end
%     results(1,i)=nbreactivation;
%     results(2,i)=breactivation;
%     rip(1,i)=size(nbspikes,2);
%     rip(2,i)=sum(nbspikes>0);
%     rip(3,i)=size(bspikes,2);
%     rip(4,i)=sum(bspikes>0);
%     
% end
% end
% 
% 
% results={};
% rip={};
% k=0;
% for i=1:size(f3.epochs{1,1},1)
%     
%     for j=1:length(f.data{1,1}{i})
%         n=k+j;
%         nbspikesind=find(f3.output{1,1}(1,n).ripple_barrierlist==0);
%     nbspikes=f3.output{1}(1,n).ripple_cell_count(nbspikesind);
%     nbreactivation=sum(nbspikes>0)/size(nbspikes,2);
%     % find all barrier
%     bspikesind=find(f3.output{1}(1,n).ripple_barrierlist==1);
%     breactivation=0;
%     bspikes=[];
%     
%     if ~isempty(bspikesind)
%         bspikes=f3.output{1}(1,n).ripple_cell_count(bspikesind);
%         breactivation=sum(bspikes>0)/size(bspikes,2);
%     end
%     
%     
%     
%     results{i}(1,j)=nbreactivation;
%     results{i}(2,j)=breactivation;
%     rip{i}(1,j)=size(nbspikes,2);
%     rip{i}(2,j)=sum(nbspikes>0);
%     rip{i}(3,j)=size(bspikes,2);
%     rip{i}(4,j)=sum(bspikes>0);
%     
%     if j==length(f3.data{1,1}{i})
%         k=k+j;
%     end
%     end
% end
%     
% 
% for i=1:size(results,2)
%     day=f3.epochs{1,1}(i,1);
%     epoch=f3.epochs{1,1}(i,2);
%     
% %     figure;
% %     nbmean=mean(results{i}(1,:));
% %     nbsem=std(results{i}(1,:))/sqrt(size(results{i}(1,:),2));
% %     bmean=mean(results{i}(2,:));
% %     bsem=std(results{i}(2,:))/sqrt(size(results{i}(2,:),2));
% %     %barweb([nbmean;bmean],[nbsem; bsem],['-barrier','+barrier'],[],[],[],[],[],[], 2, 'axis');
% %     barweb([nbmean bmean],[nbsem bsem], [], [], [], [], [], jet, [], {'-barrier';'+barrier'}, 2, 'axis');
% %     title(sprintf('1 Activation probability for %s \n day %s epoch %s',animals{1,1},num2str(day),num2str(epoch)));
%     
%     figure;
%     boxdata=[results{i}(1,:)';results{i}(2,:)'];
%     boxgroup=[zeros(size(results{i}(1,:)))'; ones(size(results{i}(2,:)))'];
%     boxplot(boxdata,boxgroup);
%     title(sprintf('2 Activation probability for %s \n day %s epoch %s',animals{1,1},num2str(day),num2str(epoch)));
%     
%     
%     
%     
% end
    
    


    
    
