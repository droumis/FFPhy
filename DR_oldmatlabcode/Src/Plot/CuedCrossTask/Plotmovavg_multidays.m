function [] = Plotmovavg_multidays(datadirectory,animalname,days,runtype,step)
% Plots the moving average for all epochs of a day
% 
% datadirectory - directory where the output from batchDIOextract is
%                 located, eg '/data14/jai/F12_/'
% animalname    -  name of the animal, eg 'F1'
% days          - vector of day or days to be analysed, eg 1
% runtype       -  'home' homebound runs 
%                  'out' outbund runs
%                  'all' all runs
% step          - number of runs for average eg 5
% 

% load file for each day

% visit =[];

% all visits across days
dayvisit=[];
% length of visit of each day across days
daylength=[0];

for i=1:size(days,2);

% open files
dayt = num2str(days(1,i));
datafilepref = datadirectory(end-5:end-1);
filename=strcat(datafilepref,num2str(days(1,i)),'.mat');
datafile = strcat(datadirectory,filename);
load(datafile);

% find no. epochs
j =  size(Data{1,days(1,i)},2);

% find data corresponding to animal

% all epochs for the day
epochmat=[];

%epochs for each day
for k = 1:j; 
    if isempty(Data{1,days(1,i)}{1,k});
        k=k+1;
    else
        runmat =[];
        if strcmp(Data{1,days(1,i)}{1,k}.Stats.animal(1,1),animalname)==1;
            % get wells
            wells=[];
            for l=1:size(Data{1,days(1,i)}{1,k}.Config,2);
                wells(1,l)=Data{1,days(1,i)}{1,k}.Config{1,l}(2,1);
            end
            % which type of run to analyse
            switch runtype;
                case 'out';
            % get all outbound runs
                runmat = Data{1,days(1,i)}{1,k}.Run(find(Data{1,days(1,i)}{1,k}.Run(:,1)==wells(1,1)),7);
                percentcorrect = Data{1,days(1,i)}{1,k}.Stats.outboundcorrect(1,1);
                case 'home';
            % get all homebound runs
                runmat = Data{1,days(1,i)}{1,k}.Run(find(Data{1,days(1,i)}{1,k}.Run(:,1)~=wells(1,1)),7);
                percentcorrect = Data{1,days(1,i)}{1,k}.Stats.inboundcorrect(1,1);
            % get all runs
                case 'all';
                runmat = Data{1,days(1,i)}{1,k}.Run(:,7);
                percentcorrect = mean(Data{1,days(1,i)}{1,k}.Run(:,7));
            end
            %epocht = num2str(k);
            % add epochs to another
            epochmat=[epochmat;runmat];
        else k=k+1;
        % add epoch to another
        end
    end   
end
   currday=size(epochmat,1);
   daylength=[daylength;currday];
   dayvisit=[dayvisit;epochmat];
 
   i=i+1;     
end

avgrunmat = moving(dayvisit,step);
            
figure;

plot(avgrunmat,'Color','red','LineWidth',4);

%title([animalname,' ','day ',dayt,' epoch ',epocht,' ', runtype, ' bound runs', char(10), 'moving average of ',num2str(step),' runs']);

%text(1,0.95,['fraction total runs correct: ',num2str(percentcorrect,2)]);

xlabel('runs');
ylabel('correct runs');



hold;
daycumsum=cumsum(daylength,1);
col=[];

o=1;
while o<size(daycumsum,1);
    currcol=[1 1 1;0.75 0.75 0.75];
    col=[col;currcol];
    o=o+2;
end
    
for n=1:size(daycumsum,1)-1;
    line([daycumsum(n,1) daycumsum(n+1,1)],[0 0],'Color',col(n,:),'LineWidth',5000)
    text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),0.05,num2str(days(1,n)),'HorizontalAlignment','center','FontSize',10);
    n=n+1;
end

plot(avgrunmat,'Color','black','LineWidth',4);

title({sprintf('Performance of %s runs for %s days %s to %s', runtype, animalname,num2str(min(days)),num2str(max(days)));...
    sprintf('%s run moving average',num2str(step))});

ylim([0 1]);
xlim([0 size(dayvisit,1)]);
line([0 max(daycumsum)],[0.5 0.5],'Color','red','LineWidth',1)
%set(gca,'XTick',0:250:max(dayvisit));


set(gca,'Layer','top');
end


