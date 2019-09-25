function [ output_args ] = PlotRerouteDist_multiday( animaldir,animalname,days )
% Plots trajectory of each run
% datadirectory - directory where the output from batchDIOextract is
% located, eg '/data14/jai/F12_/'
% animalname -  name of the animal, eg F1
% days - vector of day or days to be analysed, eg [1], [1:3]
% 


% get the current directory to go back to after the function is done
currentdir = pwd;

% get day by day data

% declare variables
% all visits across days
dayvisit=[];
% length of visit of each day across days
daylength=[0];

for d=1:size(days,2);

% open files
dayt = num2str(days(1,d));

% ---- File loading ----
% See if day number needs 0
dsz = '';
   if (d < 10)
      dsz = '0';
   end
   
% Load the mat file
datadir = '/data14/jai/'; % Specify data directory and load the file 
dayt = num2str(days(d)); % Converts the day and epoch into text
sfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'_',dsz,dayt,'.mat');
load(sfilename);
% find no. epochs
j =  size(Data{1,days(1,d)},2); 
% data for all epochs of each day
epochmat=[];

for k=1:size(Data{1,days(d)},2);
    if isempty(Data{1,days(d)}{1,k});
        k=k+1;
        
      elseif strcmp(Data{1,days(1,d)}{1,k}.Stats.epochtype,'Sleep')==1;
            k=k+1;  
        
        
    else
        if strcmp(Data{1,days(d)}{1,k}.Stats.animal(1,1),animalname)==1; %find all epochs for animal
            epoch = k;
            epocht = num2str(epoch);
            nrows=ceil(size(Data{1,days(d)}{1,epoch}.Run,1)/5);
            distance=[]; %clean matrix
            % get start and end times for each run
            for i=1:size(Data{1,days(d)}{1,epoch}.Run,1); 
                tstart=Data{1,days(d)}{1,epoch}.Run(i,3);
                tend=Data{1,days(d)}{1,epoch}.Run(i,4);
                % look up indices for time period
                [val tstartref]=min(abs((Data{1,days(d)}{1,epoch}.Pos.correcteddata(:,1)*10000-tstart)));
                [val tendref]=min(abs(Data{1,days(d)}{1,epoch}.Pos.correcteddata(:,1)*10000-tend));
                % get position values
                traj=[Data{1,days(d)}{1,epoch}.Pos.correcteddata(tstartref:tendref,2) Data{1,days(d)}{1,epoch}.Pos.correcteddata(tstartref:tendref,3) ];
                % calculate distance traveled for each run
                dist=0;
                for j=1:size(traj,1)-1;
                    dt=sqrt((traj(j+1,1)-traj(j,1))^2+(traj(j+1,2)-traj(j,2))^2);
                    dist=dist+dt;
                    j=j+1;
                end
                distance(i,1)=dist; % put values in matrix
                i=i+1;
            end
            epochmat=[epochmat;distance];
            k=k+1;
        else display(sprintf('No %s data in epoch %s ',animalname,num2str(k)));
            k=k+1;
        end
        
    end
end
    currday=size(epochmat,1);
    daylength=[daylength;currday];
    dayvisit=[dayvisit;epochmat];
 
    d=d+1;
end

  figure;



% mark each day with a different background colour


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
    text(((daycumsum(n+1,1)-daycumsum(n,1))/2+daycumsum(n,1)),9.5,num2str(days(1,n)),'HorizontalAlignment','center','FontSize',10);
    n=n+1;
end

title(sprintf('Run distance for %s days %s to %s', animalname,num2str(min(days)),num2str(max(days))));
hold on;
scatter([1:1:size(dayvisit,1)],dayvisit/100,'.','MarkerFaceColor','b','MarkerEdgeColor','b');

xlabel('Run');
ylabel('Distance (m)');
ylim([0 10]);
xlim([0 size(dayvisit,1)]);

%distlabel=[0.05*(xlimval(2)-xlimval(1))+xlimval(1) 0.90*(ylimval(2)-ylimval(1))+ylimval(1)];
%text(distlabel(1),distlabel(2),sprintf('Tot. dist. %s m',num2str(totaldist/100,3)),'FontSize',15);


end


