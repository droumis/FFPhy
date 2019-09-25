% plots the combined place fields of all place cells for a day (mean normalised rate) to see the
% distribution of placefields on the track

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

%days='[1:10]';%,'1:10';
%days = '[1:1]';
days = '[13]';




%% get placefield data
epochtype='Run';
%epochfilter = ['isequal($epochtype, ''Run'')'];
epochfilter{1} = ['isequal($epoch, 6)'];
cellfilter = '(isequal($area, ''CA1'') && ($meanrate <10 ))'  ; %excitatory cells, used runplotavgrate to see distributions for each animal
timefilter = { {'JY_getriptimes','($nripples == 0)', [], 2,'cellfilter', '(isequal($area, ''CA1''))'},{'JY_getlinvelocity', strcat('$velocity > ',num2str(minVPF))}};
f = JY_createfilter(days,'animal',animals,'days',days,'epochs',epochfilter,'cells',cellfilter,'excludetimefilter', timefilter);
%iterator = 'singlecellanal';
iterator = 'multicellanal';
f = setfilteriterator(f,iterator);
%out=setfilterfunction(f, 'JY_calcopenfieldoccupancy', {'spikes','data'});
%out=runfilter(out);
%outdata=out.output{1,1};

out2=setfilterfunction(f, 'JY_getriprate', {'ripples'});
out2=runfilter(out2);




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
    
    
    
   
    % setup figure
%     figure;
%     
%     set(gcf,'PaperUnits','inches');
%     set(gcf,'PaperSize',[2500 1000]);
%     set(gcf,'PaperPositionMode','auto');
%     set(gcf,'position',[10 10 2300 800]);
%     
    


    
    % get a list of trial times
    
    trialtimes=data{1,day}{1,epoch}.Run(:,3:4);
    
    % get a list of ripple start times
    %goodripsind=dataref(1,:);
    
%     goodrips=ripple_spikeslist{day}{epoch}(2,dataref(1,:));
%     
%     ripstart=cellfun(@(x) min(x(:,1)), goodrips)*10000;
%     
%     % find which one is not with in a run period
%     ripincl=isExcluded(ripstart,trialtimes);

barrierevents=data{1,day}{1,epoch}.Events.Barrier;

    % calculate ripples rate for each trial
    
    edges=sort([trialtimes(:,1);trialtimes(:,2)]);
    bindata=[];
    
    for i=1:size(out2.output{1,1}(1,d).nripples,2)
        if out2.output{1,1}(1,d).nripples(1,i)~=0
            bindata=[bindata;repmat(out2.output{1,1}(1,d).times(1,i),out2.output{1,1}(1,d).nripples(1,i),1)];
        end
    end
    
    %barriertrials=find(edges>
    
    ripplehist=histc(bindata.*10000,edges);
    
    trialduration=diff(edges/10000);
    
    riprate=ripplehist(1:length(trialduration),1)./trialduration;
    
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
    % gaussian smoothing
sigma=100;
width = round((6*sigma - 1)/2);
support = (-width:width);
gaussFilter = exp( -(support).^2 ./ (2*sigma^2) );
gaussFilter = gaussFilter/ sum(gaussFilter);

filtereddata=conv(out2.output{1,1}(1,d).nripples,gaussFilter,'same');

test=moving(out2.output{1,1}(1,d).nripples(1,:)',1000,'mean');


    plot(out2.output{1,1}(1,d).times*10000,out2.output{1,1}(1,d).nripples);
    plot(out2.output{1,1}(1,d).times*10000,filtereddata,'m');
    
    hold on;
    %plot(data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000, data{1,day}{1,epoch}.Pos.correcteddata(:,5),'r');
    
    % barrier label
    barrierevents=data{1,day}{1,epoch}.Events.Barrier;
    %plot
    
    plotind=1;
    

        
        if ~isempty(barrierevents);
            
                lineb=line([barrierevents(1,1) (barrierevents(1,1)+barrierevents(1,4))],[30 30],'Color','r','LineWidth',10);
 
        end

    % Print title for whole figure
    title(sprintf('Number of ripples %s for day %s epoch %s',...
        animals{1,1}, num2str(daylistmain(d,1)), num2str(daylistmain(d,2))));
    
    
%     
%         % ----Saving----
%         % Saves figure as pdf
%         % First checks if a folder called Plot exists in the processed data folder,
%         % if not, creates it.
%     
%         cd(f.animal{1,2});
%         plotdir = dir('Plot');
%         if (isempty(plotdir))
%             %an a plot folder needs to be created
%             !mkdir Plot
%         end
%     
%         % change to that directory and saves the figure with file name
%         % animal_day_epoch
%         cd(strcat(f.animal{1,2},'Plot/'));
%         figurename = strcat(animals{1,1},'ripplesvar','_',num2str(daylistmain(d,1)),'_',num2str(daylistmain(d,2)));
%     
%         saveas(gcf, figurename, 'pdf');
%     
%         % Closes the figure
%         close;
    
    
    
    
end
%clear all;









% 
% day=out.epochs{1,1}(1,1);
% 
% % collect epoch data
% daydata=[];
% epochindex=[];
% for i=1:size(out.data{1,1},2)
%     daydata=[daydata;out.data{1,1}{1,i}];
%     epochindex=[epochindex;repmat(out.epochs{1,1}(i,:),size(out.data{1,1}{1,i},1),1)];
% end
% 
% % plot combined place fields for all cells
% 
% for i=1:size(epochindex,1)
%     totalimage(:,:,i)=outdata{1,i}.normsmoothedspikerate;
% end
%     
%     totalimage(isnan(totalimage))=0;
%     imagedata=mean(totalimage,3);
%     imagesc(flipud(imagedata));
%     title(sprintf('Combined place field for all cells \n %s for day %s ',...
%         animals{1,1}, num2str(day)));
%     cd(strcat(out.animal{1,2},'Plot/'));
% figurename = sprintf('Combined_plf_%s_%s',animals{1,1},num2str(day));
% 
% 
% saveas(gcf, figurename, 'pdf');
% close;
% clear all;

