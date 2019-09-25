function [] = Plot1stvisit_multidays(datadirectory,animalname,days)
% Plots the 1st outbound visit of an animal over days
% datadirectory - directory where the output from batchDIOextract is
% located, eg '/data14/jai/F12_/'
% animalname -  name of the animal, eg F1
% days - vector of day or days to be analysed, eg [1], [1:3]
% 

% load file for each day

visit =[];
daylength=[0];
allreward=[];
for i=1:size(days,2);
currdayreward=[];
% open files
dayt = num2str(days(1,i));
datafilepref = datadirectory(end-5:end-1);
filename=strcat(datafilepref,num2str(days(1,i)),'.mat');
datafile = strcat(datadirectory,filename);
load(datafile);
   

% find no. epochs
j =  size(Data{1,days(1,i)},2);


% all epochs for the day
epochmat=[];

% find data corresponding to animal
dayvisit=[];
dayreward=[];
for k = 1:j;
    if isempty(Data{1,days(1,i)}{1,k});
        k=k+1;
    else
        outbound=[];
        reward=[];
        if strcmp(Data{1,days(1,i)}{1,k}.Stats.animal(1,1),animalname)==1;
            % get wells
            wells=[];
            
            for l=1:size(Data{1,days(1,i)}{1,k}.Config,2);
                wells(1,l)=Data{1,days(1,i)}{1,k}.Config{1,l}(2,1);
            end
           
            % get all outbound runs
            outbound = Data{1,days(1,i)}{1,k}.Run(find(Data{1,days(1,i)}{1,k}.Run(:,1)==wells(1,1)),3);
            for p=1:size(outbound,1);
                q=find(wells(1,:)==outbound(p,1));
                outbound(p,2)=q-1;
                
                p=p+1;
            end
            reward= [Data{1,days(1,i)}{1,k}.Wellinfo.rewardedwells(:,2)]';
            dayreward=[dayreward;reward];
            epochmat=[epochmat;outbound];
        else k=k+1;
        end
        
    end 
    
end
  currday=size(epochmat,1);
  currdayreward(1,1)=max(dayreward(:,1));
  currdayreward(1,2)=max(dayreward(:,2));
  allreward=[allreward;currdayreward];
  daylength=[daylength;currday];
  visit=[visit;epochmat];       
  i=i+1;        
        
end




figure;

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

cyan=[175 238 238]/255;

for n=1:size(daycumsum,1)-1;
    line([daycumsum(n,1) daycumsum(n+1,1)],[allreward(n,1) allreward(n,1)],'Color',cyan,'LineWidth',10)
    line([daycumsum(n,1) daycumsum(n+1,1)],[allreward(n,2) allreward(n,2)],'Color',cyan,'LineWidth',10)
    n=n+1;
end

hold;

plot(visit(:,2),'o',...
    'MarkerEdgeColor','r',...
    'MarkerFaceColor','r',...
    'MarkerSize',2);

xlabel('runs');
xlim([0 size(visit,1)]);
ylim([0 4]);
set(gca,'YTickLabel',{'','Well 1','Well 2', 'Well 3',''});
title({sprintf('1st visit for %s', animalname);...
    sprintf('days %s to %s',num2str(min(days)),num2str(max(days)))});
end

 
