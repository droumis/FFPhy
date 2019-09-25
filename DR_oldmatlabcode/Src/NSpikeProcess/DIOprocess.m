function [Epochdata]=DIOprocess(animaldir,day,epoch,configfile,wells,odour)
% Plots duration of DIO in line graph.
% animaldir - name of directory containing processed data
% day - day of experiment
% epoch - epoch to be analysed
% configfile - NSpike configfile for the day, stored in configdir currently
% /data14/jai/Config/
% well - vector indicating odour used eg [1 3]
% odour - odours used for the wells, specify in same order {'orange'
% 'melon'}


% get the current directory to go back to after the function is done
currentdir = pwd;

% well odour match


% ---- File loading ----
% See if day number needs 0
dsz = '';
   if (day < 10)
      dsz = '0';
   end
   
% Specify data directory and load the file 
animalname = animaldir(1:end-1);
datadir = '/data14/jai/';

% Converts the day and epoch into text
dayt = num2str(day);
epocht = num2str(epoch);
sfilename = strcat(datadir,animaldir,'/',animalname,'DIO',dsz,dayt,'.mat');

% Load the mat file
eval(['load ', sfilename]);

% Import config file for the day

configdir = '/data14/jai/Config/';

fid = fopen(strcat(configdir,configfile),'r');

allstrings=textscan(fid,'%s');
fclose(fid);
allstrings=allstrings{:};

% Extract Config from file

% Number of wells

nWellsCell = strmatch('nWells',allstrings);
nWells = str2double(allstrings(nWellsCell+1,1));

%Config = cell(nWells,1);


    findwellinfo = strmatch('well',allstrings,'exact');
    findwellinfo = [findwellinfo;size(allstrings,1)];
    for j = 1:size(findwellinfo,1)-1;
        cellstart = findwellinfo(j,1);
        cellend = findwellinfo(j+1,1)-1;
        currwell = allstrings([cellstart:cellend],1);
        well = str2double(currwell(strmatch('well',currwell,'exact')+1,1));
        inputBit = str2double(currwell(strmatch('inputBit',currwell,'exact')+1,1));
        outputBit1 = str2double(currwell(strmatch('outputBit1',currwell,'exact')+1,1));
        outputBit1Length = str2double(currwell(strmatch('outputBit1Length',currwell,'exact')+1,1));
        outputBit2 = str2double(currwell(strmatch('outputBit2',currwell,'exact')+1,1));
        outputBit2Length = str2double(currwell(strmatch('outputBit2Length',currwell,'exact')+1,1));
        outputBit3 = str2double(currwell(strmatch('outputBit3',currwell,'exact')+1,1));
        outputBit3Length = str2double(currwell(strmatch('outputBit3Length',currwell,'exact')+1,1));
        outputBit3Delay = str2double(currwell(strmatch('outputBit3Delay',currwell,'exact')+1,1));
        output0Bit = str2double(currwell(strmatch('output0Bit',currwell,'exact')+1,1));
        output0BitLength = str2double(currwell(strmatch('output0BitLength',currwell,'exact')+1,1));
    
    wellinfo=[well;inputBit; outputBit1; outputBit1Length; outputBit2; outputBit2Length; outputBit3; outputBit3Length; outputBit3Delay;output0Bit;output0BitLength];
    Config{1,j}=wellinfo;
 
    j=j+1;
    
end
 
Configinfo={'well','inputBit','outputBit1','outputBit1Length','outputBit2','outputBit2Length','outputBit3','outputBit3Length','outputBit3Delay','output0Bit','output0BitLength'}';

% ----Data extraction ----

% Specify home well

homewell = Config{1,1}(2,1);
%homewell = homewell;
% Find active wells and their next at 0

activewells =[];
for i =2:size(Config,2);
    if mean([Config{1,i}(3,1), Config{1,i}(5,1), Config{1,i}(7,1)]) ~= -1;
    curractive(1,1) = Config{1,i}(1,1);
    curractive(1,2) = Config{1,i}(2,1);
    curractive(1,3) = Config{1,i}(10,1);
    curractive(1,4) = Config{1,i}(3,1);
    activewells=[activewells;curractive];
    else activewells=activewells;
    end
    i=i+1;
    
end



Results = [];

% Get all home triggers

%Homeevents=cell(size(DIO{1,day}{1,epoch}{1,homewell}.pulsetimes(:,1),1),2);

Homeevents=[];

for i=1:size(DIO{1,day}{1,epoch}{1,homewell+1}.pulsetimes,1);
    currhome=[];
    currhome(1,1)= DIO{1,day}{1,epoch}{1,homewell+1}.pulsetimes(i,1);
    currhome(1,2)= homewell;
    Homeevents=[Homeevents;currhome];

    i=i+1;
end

% Get all home out events

homeout = [];
rewarded=activewells(find(activewells(:,3)~=-1),3);
for i=[rewarded+1]';
    
    if size(DIO{1,day}{1,epoch}{1,i}.pulsetimes,1)~=0;    
        for j=1:size(DIO{1,day}{1,epoch}{1,i}.pulsetimes,1);
            currhomeout = [DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,1),i-1];
            homeout = [homeout;currhomeout];
            j=j+1;
        end
    end
    end

    
% Match home triggers with home out events
% Out events happen within 500 time units after trigger

for i=1:size(Homeevents,1);    
    h=homeout(:,1);
    curr=Homeevents(i,1);
    k=curr-h;
    l=k(find(k<0));
    
    if size(l,1)==0;
        Homeevents(i,3)=NaN;
    elseif max(l)>=-500;
        [m,o]=find(k==max(l));
        
        Homeevents(i,3)=homeout(m,2);
    elseif max(l)<=-500
        Homeevents(i,3)=NaN;
     end
    i=i+1;
end
 

% Get all well triggers

Wellevents=[];

for i=[activewells(:,2)+1]';
    for j=1:size(DIO{1,day}{1,epoch}{1,i}.pulsetimes,1);
    currwell=[];
    currwell(1,1)= DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,1);
    currwell(1,2)= i-1;
    Wellevents=[Wellevents;currwell];
    j=j+1;
    
    end

    i=i+1;
end

% Get all well events
% Out events happen within 500 time units after trigger

for i=1:size(Wellevents,1);
    t=Wellevents(i,1);
    
    inbit=Wellevents(i,2);
    outrow=find(activewells(:,2)==inbit,1);
    outbit=activewells(outrow,4)+1;
    if outbit>0; 
    h=DIO{1,day}{1,epoch}{1,outbit}.pulsetimes(:,1);
        k=t-h;
        k=k(find(k<0));
        if size(k,1)==0;
            Wellevents(i,3)=NaN;
        elseif max(k)>=-500;
            Wellevents(i,3)=outbit-1;
        elseif max(k)<=-500
            Wellevents(i,3)=NaN;
    
        end
    else Wellevents(i,3)=NaN;
    end
        i=i+1;
end

% Pool results from home and well events
% Order according to time

Results = [Homeevents;Wellevents];

Results = sortrows(Results,1);

% Remove repeated well activations within retrigger units
list=Results(1,:);

retrigger=10000;
i=1;
while i<size(Results,1)-1 ;
    k=find(Results(:,1)>=(Results(i,1)+retrigger));
    if size(k,1)~=0;
    list=[list;Results(k(1,1),:)];
    i=k(1,1);
    else i=size(Results,1);
    end
        
end

Results=list;

% Remove repeated well activations of the same well within wellretrigger units
list=Results(1,:);

wellretrigger=50000;
i=1;
while i<size(Results,1);
    if Results(i,2)~=Results(i+1,2);
        list=[list;Results(i+1,:)];
        i=i+1;
    else    
        if  Results(i,1)+wellretrigger<Results(i+1,1);
            list=[list;Results(i+1,:)];
            i=i+1;
        else 
            i=i+1;
        end
       
    end
end

Results=list;

% Convert to text

ResultsTxt(:,1)=timetrans(Results(:,1),10000,1);

for i=1:size(Results,1);
    if Results(i,2)==homewell(1,1);
        ResultsTxt{i,2}='Home';
    else [a,b]=find(activewells(:,2)==Results(i,2));
        
        ResultsTxt{i,2}=strcat('Well',num2str(activewells(a,1)));
    
    end
    i=i+1;
end

% Rewarded or not

for i=1:size(Results,1);
    if Results(i,3)>=0;
        ResultsTxt{i,3}='Rewarded';
    else ResultsTxt{i,3}='';
    end
    i=i+1;
end


% Analyse runs

% Find first home
k=find(Results(:,2)==homewell(1,1));
kk=find(Results(k(:,1),3)>=0);
i=k(1,1);
l = find(Results(:,3)>=0);
run=[];

% Fill in form
while i<=size(Results,1)-1;

    l = find(Results(:,3)>=0);
    
    m = find(l>i);
    crun(1,1)=Results(i,2);
    crun(1,4)=Results(i,1);
    if size(m,1)~=0;
        crun(1,3)=Results(l(m(1,1),1),2);
        crun(1,5)=Results(l(m(1,1),1),1);
        crun(1,6)=crun(1,5)-crun(1,4);
        crun(1,8)=l(m(1,1),1)-i-1;
    else
        crun(1,3)=NaN;
        crun(1,5)=NaN;
        crun(1,6)=NaN;
        crun(1,8)=NaN;
    end
    if crun(1,8)==0;
        crun(1,7)=1;
    else crun(1,7)=0;
    end

    run=[run;crun];
    if size(m,1)~=0;
    i=l(m(1,1),1);
    else i=size(Results,1);
    end

end 

% Find next in sequence

for i=1:size(run,1);
    if run(i,1)==homewell(1,1);
    [a,b]=find(Results(:,1)==run(i,4));
    [c,d]=find(activewells(:,3)==Results(a,3));
    run(i,2)=activewells(c,2);
    else run(i,2)=homewell(1,1);
    end
    i=i+1;    
end


% Make text version
% Find start and end locations for each run
for i=4:6;
    Runtxt(:,i)=timetrans(run(:,i),10000,1);
end
for i=1:3;
for j=1:size(run,1);
    if run(j,i)==homewell(1,1);
        Runtxt{j,i}='Home';
    else [a,b]=find(activewells(:,2)==run(j,i));
        if size(a,1)~=0;
            Runtxt{j,i}=strcat('Well',num2str(activewells(a,1)));
        else Runtxt{j,i}='-';
        end
    j=j+1;
    end
end
i=i+i;
end




for i=1:size(run,1);
    if run(i,7)==1;
        Runtxt{i,7}='Correct';

    else
        Runtxt{i,7}='Incorrect';

    end
    
    Runtxt{i,8}=num2str(run(i,8));
    %Runtxt{i,8}=timetrans(run(i,5)/(run(i,7)+1),10000,1);
    i=i+i;
end



% Write out stats

Epochstats=[];
epoch=epoch;
duration=Results(size(Results,1),1)-Results(1,1);
outboundn=size(find(run(:,1)==homewell(1,1)),1);
inboundn=size(find(run(:,1)~=homewell(1,1)),1);
outboundcorrect=size(find(run(find(run(:,1)==homewell(1,1)),6)==1),1)/outboundn;
inboundcorrect=size(find(run(find(run(:,1)~=homewell(1,1)),6)==1),1)/inboundn;
hometrigger=size(find(Results(:,2)==homewell(1,1)),1);
homerewards=size(find(Homeevents(:,3)~=NaN),1);

name=textscan(animaldir,'%1s');
name=name{1,1};

if mod(epoch,2)==1;
    animal=strcat(name(1,1),name(2,1));
else animal=strcat(name(1,1),name(3,1));
end


odourlist=[];
for i=1:size(wells,2);
    [k l]=find(activewells(:,1)==wells(1,i));
    currod{1,1}=strcat('Well',num2str(activewells(k,1)));
    currod{1,2}=odour(1,i);
    odourlist=[odourlist;currod];
    i=1+1;
end

activewellinfo=[];

for i=1:size(activewells,1);
    activewelltriggern=size(find(Wellevents(:,2)==activewells(i,2)),1);
    currwell(1,1)=activewells(i,1);
    currwell(1,2)=activewelltriggern;
    activewellrewardsn=size(find(Wellevents(:,3)==activewells(i,4)),1);
    currwell(1,3)=activewellrewardsn;
    activewellinfo=[activewellinfo;currwell];
    i=i+1;
end

Events={};
Events.Homeevents=Homeevents;
Events.Wellevents=Wellevents;
Epochstats=[animal;epoch;duration;outboundn;inboundn;outboundcorrect;inboundcorrect;hometrigger;homerewards;{odourlist}];
Epochdata.Config=Config;
Epochdata.Configinfo=Configinfo;
Epochdata.Events=Events;
Epochdata.Stats.animal=animal;
Epochdata.Stats.epoch=epoch;
Epochdata.Stats.duration=duration;
Epochdata.Stats.outboundn=outboundn;
Epochdata.Stats.inboundn=inboundn;
Epochdata.Stats.outboundcorrect=outboundcorrect;
Epochdata.Stats.inboundcorrect=inboundcorrect;
Epochdata.Stats.hometriggers=hometrigger;
Epochdata.Stats.homerewards=homerewards;
Epochdata.Stats.odourlist=odourlist;
Epochdata.Results=Results;
Epochdata.Resultstxt=ResultsTxt;
Epochdata.Runtxt=Runtxt;
wellfinfo.activewell=activewellinfo;
wellfinfo.homewell=homewell;

Epochdata.Run=run;

fileoutdir=strcat(datadir,animaldir,'/');


% change to that directory and saves the figure with file name
% animal_day_epoch

%cd(fileoutdir);
%filename = strcat(animalname,'_',dayt,'e',epocht,'.mat');
%save(filename,'Epochdata');


%cd(currentdir);






