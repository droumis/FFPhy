function [ output_args ] = PlotRerouteVel_multiday( animaldir,animalname,daylist, barrier)
% Plots the performance for runs with different segments travelled or
% distance travelled over days
% animaldir     - directory where the output from batchDIOextract is no need 
%                 to type full path eg.'/F12_/'             
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

figure;
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
% table for collating data over days
epochstat=[];
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

epochmat=[];

for k=1:size(Data{1,daylist(1,d)},2);
    if isempty(Data{1,daylist(1,d)}{1,k});
        k=k+1;
       elseif strcmp(Data{1,daylist(1,d)}{1,k}.Stats.epochtype,'Sleep')==1;
            k=k+1;   
        
    else
        epochrunmat=[];
        if strcmp(Data{1,daylist(1,d)}{1,k}.Stats.animal(1,1),animalname)==1; %find all epochs for animal
            epoch = k;
            epocht = num2str(epoch);
            rewardedevents=Data{1,day}{1,epoch}.Events.Wellevents((Data{1,day}{1,epoch}.Events.Wellevents(:,3)~=-50),:);
            rewardedevents=rewardedevents((rewardedevents(:,4)>100),:);
            segments=[]; %clean matrix

                for i=1:size(rewardedevents,1)-1;
                    timeind=rewardedevents(i,1)+rewardedevents(i,4);
                    tstartref=min(find(Data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000>(timeind))); % for some reason +10000 seems to be accurate rather than no offset
                    timeind2=rewardedevents(i+1,1);
                    tendref=max(find(Data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000<(timeind2)));

                % get velocity info between times
                seg=abs(Data{1,daylist(1,d)}{1,epoch}.Pos.correcteddata(tstartref:tendref,5));
                segvmax=max(seg);

                %add max velocity to for each segment to table
                segments=[segments;segvmax];
                i=i+1;
                end
                % add velocities for each run to table for the epoch
                epochrunmat=[epochrunmat;segments];
       
            
            % calculate mean and sd of max velocity of each epoch
            
            
            % matrix with [day epoch mean sd]
            currepochstat=[];
            currepochstat(1,1)=daylist(1,d);
            currepochstat(1,2)=epoch;
            currepochstat(1,3)=mean(epochrunmat);
            currepochstat(1,4)=std(epochrunmat);
            epochstat=[epochstat;currepochstat];
            
            k=k+1;
        else display(sprintf('No %s data from day %s epoch %s ',animalname,num2str(daylist(1,d)),num2str(k)));
            k=k+1;
        end
    
        
    end
    
end


errorbar(epochstat(:,3),epochstat(:,4),'.k');
ylim([0 60]);
xlabel('Run Epoch');

title({sprintf('Maximum velocity for %s days %s to %s',animalname,num2str(min(daylist)),num2str(max(daylist)))});




end


