function [ disttimetable ] = PlotRerouteSeg_multiday( animaldir,animalname,daylist, barrier)
% Plots the duration for runs with different segments travelled or
% distance travelled over days
% animaldir     - directory where the output from batchDIOextract is no need 
%                 to type full path eg.'F12_'             
% animalname    - name of the animal, eg 'F1'
% daylist       - vector of day or days to be analysed, eg [1], [1:3]
% barrier       - vector containing epoch and path on which barriers are placed
%                 first number indicates the day,second indicates the epoch for that day 
%                 and third specify the intersections between which a
%                 barrier is placed
%                 eg. [2 2 1 2]
% needs connectivity matrix located in
% /home/jai/Src/NSpikeProcess/rerouteconnectivity.mat
% assumes partial segment crossings as 40cm (see segment identity 9 of distance table in connectivity matrix)

% set moving window for running average
movwin=10;
% get the current directory to go back to after the function is done
currentdir = pwd;
% get day by day data
% declare variables
% all visits across days
dayvisit=[];
% length of visit of each day across days
daylength=[0 0];
% stats for each epoch of each day [day epoch run_number]
dayepochstat=[];
% load connectivity matrix and distance table
load('/home/jai/Src/NSpikeProcess/rerouteconnectivity.mat');
rcon=rerouteconnectivity;
clear rerouteconnectivity;
for d=1:size(daylist,2);
% open files
dayt = num2str(daylist(1,d));
% ---- File loading ----
% See if day number needs 0
dsz = '';
   if (d < 10)
      dsz = '0';
   end
% Load the mat file
datadir = '/data14/jai/'; % Specify data directory and load the file 
dayt = num2str(daylist(1,d)); % Converts the day and epoch into text
sfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'_',dsz,dayt,'.mat');
load(sfilename);
posfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'_',dsz,dayt,'linpos.mat');
load(posfilename);
% find no. epochs
j =  size(Data{1,daylist(1,d)},2); 
% data for all epochs of each day
epochrunmat=[];
epochmat=[];
epochstat=[];
for k=1:size(Data{1,daylist(1,d)},2);
    if isempty(Data{1,daylist(1,d)}{1,k});
        k=k+1;
    else
        if strcmp(Data{1,daylist(1,d)}{1,k}.Stats.animal(1,1),animalname)==1; %find all epochs for animal
            epoch = k;
            epocht = num2str(epoch);
            rewardedevents=Data{1,day}{1,epoch}.Events.Wellevents((Data{1,day}{1,epoch}.Events.Wellevents(:,3)~=-50),:);
            rewardedevents=rewardedevents((rewardedevents(:,4)>100),:);
            segments=[]; %clean matrix
            % get start and end times for each run
%             for i=1:size(Data{1,daylist(1,d)}{1,epoch}.Run,1); 
%                 tstart=Data{1,daylist(1,d)}{1,epoch}.Run(i,3);
%                 tend=Data{1,daylist(1,d)}{1,epoch}.Run(i,4);
%  
                for i=1:size(rewardedevents,1)-1;
                    timeind=rewardedevents(i,1)+rewardedevents(i,4);
                    tstartref=min(find(Data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000>(timeind))); % for some reason +10000 seems to be accurate rather than no offset
                    timeind2=rewardedevents(i+1,1);
                    tendref=max(find(Data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000<(timeind2)));


                % get segment info from linpos file
%                 startind=find(linpos{1,day}{1,epoch}.statematrix.time(:,1)*10000>tstart);
%                 tstartref=min(startind);
%                 endind=find(linpos{1,day}{1,epoch}.statematrix.time(:,1)*10000<tend);
%                 tendref=max(endind);
                % get segment info between times
                seg=linpos{1,day}{1,epoch}.statematrix.segmentIndex(tstartref:tendref,1);
                % get exit times during period
                
                estartref=find(linpos{1,day}{1,epoch}.trajmatrix(:,1)*10000>timeind+10000); %removed +10000
                if ~isempty(estartref);
                estartref=estartref(1,1);
                [val eendref]=min(abs(linpos{1,day}{1,epoch}.trajmatrix(:,1)*10000-timeind2));
                if ~isempty(estartref);
                exitref=linpos{1,day}{1,epoch}.trajmatrix(estartref(1,1):eendref,:);
                else i=size(Data{1,daylist(1,d)}{1,epoch}.Run,1); 
                end
                % concatenate columns to generate integer representation of
                % path
                for q=1:size(exitref,1);
                    exitref(q,4)=str2num(sprintf('%d%d',exitref(q,2),exitref(q,3)));
                    
                end
                % remove segments containing barrier
                for b=1:size(barrier,1);
                    if barrier(b,1)==daylist(1,d);
                        if barrier(b,2)==epoch;
                            fdirect=str2num(sprintf('%d%d',barrier(b,3),barrier(b,3)));
                            forward=find(exitref(:,4)~=fdirect);
                            exitref=exitref(forward,:);
                            rdirect=str2num(sprintf('%d%d',barrier(b,4),barrier(b,4)));
                            reverse=find(exitref(:,4)~=rdirect);   
                            exitref=exitref(reverse,:);
                        end
                    end
                    b=b+1;
                end
                
                % find segment number and distance from connectivity matrix rcon
                for q=1:size(exitref,1);
                    exitref(q,5)=rcon(exitref(q,2),exitref(q,3));
                    exitref(q,6)=disttb(exitref(q,5),1);
                end
                % remove incomplete path crossings
                exitref=exitref(find(exitref(:,2)-exitref(:,3)~=0),:);
                segments(i,1)=size(exitref,1); % put values in matrix
                segments(i,2)=sum(exitref(:,6)); % distance for path
                segments(i,3)=(linpos{1,day}{1,epoch}.statematrix.time(tendref,1)-linpos{1,day}{1,epoch}.statematrix.time(tstartref,1))*10000;
                segments(i,4)=max(Data{1,day}{1,epoch}.Pos.correcteddata(tstartref:tendref,5));
                
                if segments(i,1)==0;
                    display(sprintf('Respecify well location for day %d epoch %d reward index %d start at %d end at %d',day,epoch,i,tstartref,tendref));
                end
                i=i+1;
                end
                end
            epochrunmat=[epochrunmat;size(Data{1,daylist(1,d)}{1,epoch}.Run,1)];
            epochmat=[epochmat;segments];
            % matrix with [day epoch run_number]
            currepochstat=[];
            currepochstat(1,1)=daylist(1,d);
            currepochstat(1,2)=epoch;
            currepochstat(1,3)=size(Data{1,daylist(1,d)}{1,epoch}.Run,1);
            epochstat=[epochstat;currepochstat];
            
            k=k+1;
        else display(sprintf('No %s data from day %s epoch %s ',animalname,num2str(daylist(1,d)),num2str(k)));
            k=k+1;
        end
    
        
    end
    
end
    %currday=size(epochmat,1);
    if isempty(epochrunmat);
        d=d+1;
    else
    currday(1,1)=sum(epochrunmat,1); % number of runs for the day
    currday(1,2)=sum(epochmat(:,2)); % distance for the day
    daylength=[daylength;currday]; % number of runs for each day
    dayvisit=[dayvisit;epochmat]; % number of segments for each run over days
    dayepochstat=[dayepochstat;epochstat];
    d=d+1;
    end
end


% mark each day with a different background colour
daycumsum=cumsum(daylength,1);
col=[];
o=1;
while o<size(daycumsum,1);
    currcol=[1 1 1;0.75 0.75 0.75];
    col=[col;currcol];
    o=o+2;
end


figure;

% generate a table listing all the distances and average time taken

distancelist=unique(dayvisit(:,2));
disttimetable=[]; %table [pathdistance meantime sdtime n maxvel sdmaxvel]

for d=1:size(distancelist);
    currdtt(1,1)=distancelist(d,1);
    currdtt(1,2)=mean(dayvisit(find(dayvisit(:,2)==distancelist(d,1)),3));
    currdtt(1,3)=std(dayvisit(find(dayvisit(:,2)==distancelist(d,1)),3));
    currdtt(1,4)=size(find(dayvisit(:,2)==distancelist(d,1)),1);
    currdtt(1,5)=mean(dayvisit(find(dayvisit(:,2)==distancelist(d,1)),4));
    currdtt(1,6)=std(dayvisit(find(dayvisit(:,2)==distancelist(d,1)),4));
    disttimetable=[disttimetable;currdtt];
    d=d+1;
end

hist(dayvisit(dayvisit(:,2)==222,3),100000);

% plot duration for runs of difference distances

figure;

disttimetable=disttimetable(disttimetable(:,4)>1,:);

errorbar(disttimetable(:,1)/100,disttimetable(:,2)/10000,disttimetable(:,3)/10000,'xk');

xlabel('distance of path (m)');
ylabel('time (s)');
title(sprintf('%s travel time for day %s to %s',animalname,num2str(min(daylist)),num2str(max(daylist))));

% plot max speed

figure;

disttimetable=disttimetable(disttimetable(:,4)>1,:);

errorbar(disttimetable(:,1)/100,disttimetable(:,5),disttimetable(:,6),'xk');

xlabel('distance of path (m)');
ylabel('velocity (cm/s)');
title(sprintf('%s max velocity for day %s to %s',animalname,num2str(min(daylist)),num2str(max(daylist))));

end


