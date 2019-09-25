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

days='[10]';%,'1:10';



%epochfilter{1} = ['($switchday > 0) & ($tasknum == 1)']; %just analyze days where switching between tasks
%epochfilter{2} = ['($switchday > 0) & ($tasknum == 2)'];
%epochfPF = ['($switchday > 0)'];

epochfilter{1} = ['isequal($epochtype, ''Run'')'];
%epochfilter{1} = ['isequal($epoch, 2)'];
cellfilter = '(isequal($area, ''CA1'') && ($meanrate <7))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
%timefilter = { {'JY_getlinvelocity', '$velocity <2'} };
timefilter = { {'JY_getriptimes','($nripples ==0)', [], 3,'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}, {'JY_getlinvelocity', '$velocity <0.01'}};
%timefilter2 = { {'getriptimes',[], [],'cellfilter', '(isequal($area, ''CA1''))'}};
%timefilter2 = { {'getriptimes','($nripples > 0)',[],'minthresh',2,'cellfilter',cellfilter,'tetfilter',[1 2 3 4]}};

%timefilter = { {'JY_getlinvelocity', '(($velocity) >= 0))', 6} };
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);



%only include cells with placefields
%if minPeakPF>0
%    includecells = calcincludecells(minVPF, minPeakPF, animals, epochfPF, lessthan);
%5    f = excludecellsf(f, includecells);
%end
%-----------------------------------------------------------


%run function- single cells
%--------------------------------------------
iterator = 'singlecellanal';

f = setfilteriterator(f,iterator);

out=setfilterfunction(f, 'JY_calcopenfieldoccupancy', {'spikes','data'});

out=runfilter(out);

% place fields are stored here
outdata=out.output{1,1};


% Plotting

daylistmain=unique(f.epochs{1,1}(:,1));
animaldir='I1_';
datadir = '/data14/jai/';
animalname = animals{1,1};
daylistmain=unique(f.epochs{1,1}(:,[1 2]),'rows');

for d=1:size(daylistmain,1);
    
    day=daylistmain(d,1)
    epoch=daylistmain(d,2);
    
    dsz = '';
    if (day < 10)
        dsz = '0';
    end
    
    dayt = num2str(day);
    
    % load place field file (plf)
    %     filename=strcat(datadir,animaldir,'/',animalname,'plf',dsz,dayt,'.mat');
    %     load(filename);
    
    % Load the mat file, needed for barrier information
    sfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'data',dsz,dayt,'.mat');
    load(sfilename);
    
    posfilename = strcat(datadir,animaldir,'/',animals{1,1},'linpos',dsz,dayt,'.mat');
    load(posfilename);
    
    % get rewarded wells
    Wells=data{1,day}{1,epoch}.Wellinfo.rewardedwells;
    % get reward well positions
    Wellpos=linpos{1,day}{1,epoch}.wellSegmentInfo.wellCoord;
    
    
    % barrier label
    barrierevents=data{1,day}{1,epoch}.Events.Barrier;
    
    
    
    % load spike summary file
    filename=strcat(f.animal{2},f.animal{1},'ripple_cell_count',dsz,num2str(day),'.mat');
    load(filename);
    % count of active cells in each ripple
    ripplecellcount=ripple_cell_count{day}{epoch};
    % find all ripples with more than one cell active
    rippleind=find(ripplecellcount>2);
    
    % setup figure
%     figure;
%     
%     set(gcf,'PaperUnits','inches');
%     set(gcf,'PaperSize',[2500 1000]);
%     set(gcf,'PaperPositionMode','auto');
%     set(gcf,'position',[10 10 2300 800]);
%     
    
    % find how many plots are needed
    tmpdata={};
    rippcounter=0;
    dataref=[];
    for plotind=1:length(rippleind);
        
        % get a list of cells active in each ripple
        % spikes are found in ripple_spikeslist
        uniquecelllist=unique(ripple_spikeslist{day}{epoch}{2,rippleind(plotind)}(:,[8 9]),'rows');
        cellsplacefield={};
        cellcounter=0;
        
        % get the corresponding place field of each cell
        
        for cellind=1:size(uniquecelllist,1);
            
            % get place field information
            %             if max(celldata{1,day}{1,epoch}{1,uniquecelllist(cellind,1)}{1,uniquecelllist(cellind,2)}.plf(:))<100;
            
            % find epoch
            % find tet
            currtetind=uniquecelllist(cellind,1);
            % find cell
            currcellind=uniquecelllist(cellind,2);
            % find corresponding epoch index in f.data
            epochind=find(f.epochs{1,1}(:,2)==epoch);
            % find tet/cell reference in f.data
            currdataind=find(f.data{1,1}{1,epochind}(:,1)==currtetind & f.data{1,1}{1,epochind}(:,2)==currcellind);
            if ~isempty(currdataind)
                
                % find corresponding placefield in outdata
                outdataind=(epochind-1)*size(f.data{1,1}{1,epochind},1)+currdataind;
                
                cellsplacefield{1,cellind}=outdata(1,outdataind).smoothedspikerate;
                binx=outdata(1,outdataind).xticks;
                biny=outdata(1,outdataind).yticks;
                %cellsplacefield{1,cellind}(isnan(cellsplacefield{1,cellind}))=0;
                
                cellcounter=cellcounter+1;
            else
                cellind=cellind+1;
            end
            
            %else cellsplacefield{1,cellind}=zeros(size(celldata{1,day}{1,epoch}{1,uniquecelllist(cellind,1)}{1,uniquecelllist(cellind,2)}.plf_norm));
            
            %end
            %cellind=cellind+1;
            
        end
        
        imagedata=max(cat(3,cellsplacefield{:}),[],3);
        imagedatamin=min(cat(3,cellsplacefield{:}),[],3);
        imagedata(imagedatamin==-1)=-1;
        
        if ~sum(imagedata(:))==0;
            tmpdata{rippcounter+1}=imagedata; % store current image in tmpdata
            dataref(1,rippcounter+1)=rippleind(plotind);
            dataref(2,rippcounter+1)=cellcounter;
            rippcounter=rippcounter+1;
            plotind=plotind+1;
            
        else
            plotind=plotind+1;
            
        end
    end
    
    % get 10 ripples at a time and compute variance for each pixel across
    % the 10 ripples
    movwin=11;
    index=[];
    index(:,1)=[1:1:(length(tmpdata))]';
    index(:,2)=index(:,1)+movwin-1;
    index(index(:,2)>(length(tmpdata)-1),2)=length(tmpdata);
    
    imagedata={};
    imagevar=[];
    % calculate variance of variance of 10 trial window
    for i=1:size(index,1);
        tempmat =tmpdata(index(i,1):index(i,2));
        imagedata{i}=var(cat(3,tempmat{:}),[],3);
        imagevar(i)=var(imagedata{i}(:));
    end
    
    imagevar=imagevar./max(imagevar);
    
    
    % compare using 2d crosscor euclidean distance between two images
    
    dist=[];
    
    for i=1:(length(tmpdata)-1);
        tempdist=pdist2(tmpdata{i},tmpdata{i+1});
        dist(i)=mean(tempdist(:));
    end
    
    % generate random draws to find chance
    
%     randind=randi(length(tmpdata),[100 length(tmpdata)]);
%     for i=1:size(randind,1);
%         for j=1:(size(randind,2)-1);
%            tempdist=pdist2(tmpdata{randind(i,j)},tmpdata{randind(i,j+1)});
%         testdist(i,j)=mean(tempdist(:));
%         end
%         testdist(i,:)=testdist(i,:)./max(testdist(i,:));
%     end 
%     
%     meantestdist=mean(testdist(1:74,:),1);
%     stdtestdist=std(testdist(1:74,:),1,1);
%     setestdist=stdtestdist./sqrt(74);
%         
    dist=dist./max(dist);
  
    
    % get a list of trial times
    
    trialtimes=data{1,day}{1,epoch}.Run(:,3:4);
    
    % get a list of ripple start times
    %goodripsind=dataref(1,:);
    
    goodrips=ripple_spikeslist{day}{epoch}(2,dataref(1,:));
    
    ripstart=cellfun(@(x) min(x(:,1)), goodrips)*10000;
    
    % find which one is not with in a run period
    ripincl=isExcluded(ripstart,trialtimes);
    
    
    figure;
    
    % plot trials
    
    linen=line([min(trialtimes(:,1)) max(trialtimes(:,2))],[0 0],'Color','green','LineWidth',5000);
    hold on;
    col=[];
    o=1;
    while o<size(trialtimes,1);
        currcol=[1 1 1;0.75 0.75 0.75];
        col=[col;currcol];
        o=o+2;
    end
    
    for n=1:size(trialtimes,1)-1;
        linen=line([trialtimes(n,1) trialtimes(n,2)],[0 0],'Color',col(n,:),'LineWidth',5000);
        hold on;
    end
    
    %plot(ripstart,imagevar,'-bx');
    hold on;
    plot(ripstart(1:end-1),dist,'-kx','LineWidth',1);
    
    hold on;
    
    % barrier label
    barrierevents=data{1,day}{1,epoch}.Events.Barrier;
    %plot
    
    plotind=1;
    
    for plotind=1:length(tmpdata);
        
        if ~isempty(barrierevents);
            
            % get run start and run end times
            TS=min(ripple_spikeslist{day}{epoch}{2,dataref(1,plotind)}(:,1))*10000;
            TE=max(ripple_spikeslist{day}{epoch}{2,dataref(1,plotind)}(:,1))*10000;
            BS=barrierevents(:,1);
            BE=barrierevents(:,1)+barrierevents(:,4);
            barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
            if sum(barriertest)>0;
                plot(ripstart(plotind),0.5,'.r');
                
                
            end
            
        end
        
        
    end
    % Print title for whole figure
    title(sprintf('Variance of combined place fields during ripples of %s for day %s epoch %s \n moving window of %s ripples',...
        animals{1,1}, num2str(daylistmain(d,1)), num2str(daylistmain(d,2)),num2str(movwin)));
    
    
    
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
        figurename = strcat(animals{1,1},'ripplesvar','_',num2str(daylistmain(d,1)),'_',num2str(daylistmain(d,2)));
    
        saveas(gcf, figurename, 'pdf');
    
        % Closes the figure
        close;
    
    
    
    
end
clear all;



























