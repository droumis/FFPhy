Veqn = '>=0'
minV =  str2num(Veqn(end))
maxstage = 3% [1 2 3]
minVPF = 2 %cm/sec
minPeakPF = 3
lessthan=0
includestates = 6

%Animal selection
%-----------------------------------------------------
animals = {'I1'};
%animals = {'Barack', 'Calvin', 'Dwight'};
%-----------------------------------------------------




%Filter creation
%--------------------------------------------------------
% day filter

days='[6]';%,'1:10';



%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

epochfilter{1} = ['isequal($epochtype, ''Run'')'];




cellfilter = '(isequal($area, ''CA1'') && ($meanrate >0))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '$velocity <2'} };
timefilter = { {'JY_getriptimes','($nripples ==0)', [], 3,'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { {'JY_getlinvelocity', '(($velocity) >= 0))', 6} };
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);

% pool all spike data for each epoch

animaldir='I1_';
datadir = '/data14/jai/';
for i=1:length(f.epochs{1,1});
    % epoch defined in f.epoch
    day=f.epochs{1,1}(i,1);
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    
    epoch=f.epochs{1,1}(i,2);
    % cells defined in f.data
    epochspikes=[];
    for j=1:size(f.data{1,1}{1,i},1);
        tet=f.data{1,1}{1,i}(j,1);
        cell=f.data{1,1}{1,i}(j,2);
        
        % get cell data
        spkfilename = strcat(datadir,animaldir,'/',animals,'spikes',dsz,num2str(day),'.mat');
        eval(['load ', spkfilename{1,1}]);
        currspikes=spikes{1,day}{1,epoch}{1,tet}{1,cell}.data;
        currspikes(:,8)=tet;
        currspikes(:,9)=cell;
        epochspikes=[epochspikes;currspikes];
    end
    
    % get included times

if ~isempty(excludetime)
    exclind = isExcluded(epochspikes(:,1),excludetime);  %find excluded time indices
    %pospos = posposorig(exclind==1,:);
    % get exclude index
    
    goodspikes = epochspikes(exclind==0,:);
    badspikes=epochspikes(exclind==1,:);

    
end
    epochspikes=goodspikes;
    
    uniquecelllist=unique(epochspikes(:,[8 9]),'rows');
    
    %Plot time vs spikes
    
    % first plot exlude time period (equals ripples)
    % excluded time defined in f.excluded times
    excludetime=f.excludetime{1,1}{1,i};
    
    for inds=1:size(excludetime,1)
        line([excludetime(inds,1) excludetime(inds,2)],[0 0],'Color',[0.5 0.5 0.5],'LineWidth',5000);
        hold on;
    end
    
    cmap=hsv(size(uniquecelllist,1));
    
    for cellind=1:size(uniquecelllist,1);
        
        
        currcellspikes=epochspikes(epochspikes(:,8)==uniquecelllist(cellind,1) & epochspikes(:,9)==uniquecelllist(cellind,2),:);
        

        plot(currcellspikes(:,1),cellind,'.', 'MarkerFaceColor',cmap(cellind,:),'MarkerEdgeColor',cmap(cellind,:));
        hold on;
    end
        
        ylim([0 size(uniquecelllist,1)+1]);
        
        % excluded time defined in f.excluded times
        excludetime=f.excludetime{1,1}{1,i};

        PlotopenfieldrateReroute(animaldir,tet,cell,day,epoch, 2,1,0,0,1,0,excludetime);
    end