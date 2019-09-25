function [Epochdata]=RerouteProcess(animaldir,animalnames,day,epoch,configfile,pixdim,barrierport,runepoch,sleepepoch,ranges,rewardedwells)
% Generates data for reroute task
% includes lens distortion correction
% animaldir - name of directory containing processed data
% animalnames - array with animals {'H1'} for one animal or {'H1'; 'H2'} for
% alternating runs
% day - day of experiment
% epoch - epoch to be analysed
% configfile - NSpike configfile for the day, stored in configdir currently
% /data14/jai/Config/
% --------------------------------
% pixdim - cm per pixel see batchDIOextractReroute
% --------------------------------
% barrierport - output port number used to braise barrier
% runepoch -  vector of run epochs
% sleepepoch - vector of sleep epochs
% ranges - matrix of times for each epoch from ranges in times.mat


%% file loading
% get the current directory to go back to after the function is done
currentdir = pwd;

if (nargin < 6);
    barrierport=[];
end

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


% name of animal for each epoch

if size(animalnames,1)==1; % if only one name is provided, name all epochs the same
    epochanimal=animalnames;
else
    if mod(epoch,2)==1;
        epochanimal=animalnames(1); % otherwise, all odd epochs will have the first name
    else epochanimal=animalnames(2); % all even epochs will have the second name
    end
    
end

% load times.mat
times=load(strcat(datadir,animalname,'/',animalname,dsz,dayt,'/times.mat'));



%% Extract DIO settings from Config file
% DIO information is only processed for run epochs
% sleep and run epochs are defined in times.mat and by
% batchDIOextractReroute.m

% Number of wells

nWellsCell = strmatch('nWells',allstrings);
nWells = str2double(allstrings(nWellsCell+1,1));

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

if any(runepoch==epoch);
    epochtype='Run';
    
    
    
    %% Find when barriers are raised
    % find port raising barrier, port specified by barrierport
    barrierevents=[];
    if ~isempty(barrierport);
        bpn=barrierport+1; % offsets matlab matrix reference
        j=1;
        
        while j<=size(DIO{1,day}{1,epoch}{1,bpn}.pulsetimes,1);
            if j==size(DIO{1,day}{1,epoch}{1,bpn}.pulsetimes,1);
                currwell=[];
                currwell(1,1)= DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(j,1);
                currwell(1,2)= bpn-1;
                currwell(1,4)=DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(j,2)-DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(j,1);
                barrierevents=[barrierevents;currwell];
                j=j+2;
            else
                
                currwell=[];
                currwell(1,1)= DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(j,1);
                currwell(1,2)= bpn-1;
                currwell(1,4)=DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(j,2)-DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(j,1);
                k=j;
                while k<=size(DIO{1,day}{1,epoch}{1,bpn}.pulsetimes,1)-1;
                    if DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(k,2)==DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(k+1,1);
                        if k==size(DIO{1,day}{1,epoch}{1,bpn}.pulsetimes,1)-1;
                            currwell(1,4)=DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(k,2)-DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(j,1);
                            barrierevents=[barrierevents;currwell];
                            k=k+1;
                            j=k+1;
                        else
                            k=k+1;
                        end
                    else
                        currwell(1,4)=DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(k,2)-DIO{1,day}{1,epoch}{1,bpn}.pulsetimes(j,1);
                        barrierevents=[barrierevents;currwell];
                        j=k+1;
                        k=size(DIO{1,day}{1,epoch}{1,bpn}.pulsetimes,1);
                    end
                end
            end
            
        end
    end
    
    
    
    %% Find active wells
    
    activewells =[];
    for i =1:size(Config,2);
        if mean([Config{1,i}(3,1), Config{1,i}(5,1), Config{1,i}(7,1)]) ~= -1;
            curractive(1,1) = Config{1,i}(1,1); % well number -1
            curractive(1,2) = Config{1,i}(2,1); % well input bit -1
            curractive(1,3) = Config{1,i}(3,1); % well out put bit
            activewells=[activewells;curractive];
        else activewells=activewells;
        end
        i=i+1;
        
    end
    
    Results = [];
    
    %% Get all reward events from DIO output bits
    
    % Get all well events
    % Out events happen within 500 time units after trigger
    % currently specific output bit time triggers assumed to occur after input
    % time triggers and are not explicitly stated in well events
    % remove repeated trigger artifacts
    % if the end time of current output is equal to the start time of the
    % next output, they belong to the same output
    
    outputtriggertimes=[];
    retriggertime=5;
    
    for i=[activewells(:,3)+1]';
        
        % make matrix with 1st column = start times; second column =end
        % times but shifted one
        shiftt=[];
        shiftt(:,1)=[DIO{1,day}{1,epoch}{1,i}.pulsetimes(:,1);0];
        shiftt(:,2)=[0;DIO{1,day}{1,epoch}{1,i}.pulsetimes(:,2)];
        shiftt(:,3)=abs(shiftt(:,1)-shiftt(:,2)); % compares the start and end of each pulse with the next one
        
        % get all pulses that are too close to eachother as defined by retriggertime or start
        % immediately after the previous one
        
        badpulseind=find(shiftt(:,3)<=retriggertime);
        j=1;
        while  j<=size(DIO{1,day}{1,epoch}{1,i}.pulsetimes,1);
            if ~ismember(j+1,badpulseind);
                currwell=[];
                currwell(1,1)= DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,1);
                currwell(1,2)= i-1;
                currwell(1,4)=DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,2)-DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,1);
                outputtriggertimes=[outputtriggertimes;currwell];
                j=j+1;
            else
                % find the next pulse after the bad pulse that is after the
                % trigger time
                nextgoodpulseind=find(shiftt(:,3)>retriggertime);
                nextgoodpulseind=nextgoodpulseind(nextgoodpulseind>j);
                nextgoodpulse=nextgoodpulseind(1)-1;
                currwell=[];
                currwell(1,1)= DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,1);
                currwell(1,2)= i-1;
                currwell(1,4)=DIO{1,day}{1,epoch}{1,i}.pulsetimes(nextgoodpulse,2)-DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,1);
                outputtriggertimes=[outputtriggertimes;currwell];
                j=nextgoodpulse+1;
            end
            
        end
    end
    
    
    
    if ~isempty(outputtriggertimes);
        
        outputtriggertimes = sortrows(outputtriggertimes,1);
    end
    
    
    
    %% get all well activations from input DIO bits
    
    Wellevents=[];
    retriggertime=1000;
    inputtriggertimes=[];
    %last_combined_flag = 0;
    
    for i=[activewells(:,2)+1]';
        
        % make matrix with 1st column = start times; second column =end
        % times but shifted one
        shiftt=[];
        shiftt(:,1)=[DIO{1,day}{1,epoch}{1,i}.pulsetimes(:,1);0];
        shiftt(:,2)=[0;DIO{1,day}{1,epoch}{1,i}.pulsetimes(:,2)];
        shiftt(:,3)=abs(shiftt(:,1)-shiftt(:,2)); % compares the start and end of each pulse with the next one
        
        % get all pulses that are too close to eachother as defined by retriggertime or start
        % immediately after the previous one
        
        badpulseind=find(shiftt(:,3)<=retriggertime);
        j=1;
        while  j<=size(DIO{1,day}{1,epoch}{1,i}.pulsetimes,1);
            if ~ismember(j+1,badpulseind);
                currwell=[];
                currwell(1,1)= DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,1);
                currwell(1,2)= i-1;
                currwell(1,4)=DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,2)-DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,1);
                inputtriggertimes=[inputtriggertimes;currwell];
                j=j+1;
            else
                % find the next pulse after the bad pulse that is after the
                % trigger time
                nextgoodpulseind=find(shiftt(:,3)>retriggertime);
                nextgoodpulseind=nextgoodpulseind(nextgoodpulseind>j);
                nextgoodpulse=nextgoodpulseind(1)-1;
                currwell=[];
                currwell(1,1)= DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,1);
                currwell(1,2)= i-1;
                currwell(1,4)=DIO{1,day}{1,epoch}{1,i}.pulsetimes(nextgoodpulse,2)-DIO{1,day}{1,epoch}{1,i}.pulsetimes(j,1);
                inputtriggertimes=[inputtriggertimes;currwell];
                j=nextgoodpulse+1;
            end
            
        end
    end
    
    
    if ~isempty(inputtriggertimes);
        
        inputtriggertimes = sortrows(inputtriggertimes,1);
    end
    
    
    
    % Get all well events
    % Out events happen within 500 time units after trigger
    % currently specific output bit time triggers assumed to occur after input
    % time triggers and are not explicitly stated in well events
    
    %% Combine input and output into one big matrix
    Wellevents=inputtriggertimes;
    
    for i=1:size(Wellevents,1);
        t=Wellevents(i,1);
        
        inbit=Wellevents(i,2);
        outrow=find(activewells(:,2)==inbit,1);
        outbit=activewells(outrow,3)+1;
        if outbit>0;
            if ~isempty(outputtriggertimes);
                h=outputtriggertimes(:,1);
                k=t-h;
                k=k(find(k<0));
                if size(k,1)==0;
                    Wellevents(i,3)=-50;
                elseif max(k)>=-500;
                    Wellevents(i,3)=outbit-1;
                elseif max(k)<=-500
                    Wellevents(i,3)=-50;
                    
                end
            else Wellevents(i,3)=-50;
            end
        else Wellevents(i,3)=-50;
        end
        i=i+1;
    end
    
    Wellevents=sortrows(Wellevents,1);
    
    % Pool results from home and well events
    % Order according to time
    
    Results = Wellevents;
    Results(:,5)=Results(:,1)+Results(:,4);
    
    Results = sortrows(Results,1);
    
    % Convert to text
    
    ResultsTxt(:,1)=timetrans(Results(:,1),10000,1);
    ResultsTxt(:,4)=timetrans(Results(:,4),10000,1);
    ResultsTxt(:,5)=timetrans(Results(:,5),10000,1);
    
    
    for i=1:size(Results,1);
        [a,b]=find(activewells(:,2)==Results(i,2));
        
        ResultsTxt{i,2}=strcat('Well',num2str(activewells(a,1)));
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
    
    %% Analyse runs to generate a summary of events in epoch
    
    % find start of epoch, defined by the start time
    
    epochstart=times.ranges(epoch+1,1);
    epochend=times.ranges(epoch+1,2);
    
    run=[];
    
    
    if mean(Results(:,3))~=-50 % only process when a reward is triggered
        
        
        % first trail is defined as the time between the start of the epoch and
        % the first correct reward well
        % since the the wells are preloaded, sometimes the 1st well trigger
        % does not release additional reward.
        
        % find out the first trigger
        % compare first recorded well trigger with reward delivery, if the
        % first well trigger is the same, then ignore, if not, register
        
        % first reward deliver index
        
        
        firstreward=find(Results(:,3)~=-50,1,'first');
        
        if length(unique(Results(1:firstreward,2)))~=1
            firstindex=find(Results(:,4)>1000,1,'first');
            secondindex=firstreward;
            
            
            
            
            % first trial
            crun=[];
            crun(1,1)=-33; % trial start well
            crun(1,3)=epochstart; % trial start time
            crun(1,2)=Results(firstindex,2); % trial end well
            crun(1,4)=Results(firstindex,1); % trial end time
            crun(1,5)=crun(1,4)-crun(1,3);
            
            %2nd trial
            
            crun(2,1)=Results(firstindex,2); % trial start well
            crun(2,3)=Results(firstindex,5); % trial start time
            crun(2,2)=Results(secondindex,2); % trial end well
            crun(2,4)=Results(secondindex,1); % trial end time
            crun(2,5)=crun(2,4)-crun(2,3);
            
            
            % test if run contains a barrier
            % looks for overlap between run time and barrier times
            % get run start and run end times
            
            for ii=size(crun,1)
                
                TS=crun(ii,3);
                TE=crun(ii,4);
                
                % get barrier on times
                if ~isempty(barrierevents);
                    BS=barrierevents(:,1);
                    BE=barrierevents(:,1)+barrierevents(:,4);
                    barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
                    crun(ii,6)=sum(barriertest)>0;
                else crun(ii,6)=0;
                end
            end
            
            run=crun;
            
        else
            % first trial
            crun=[];
            crun(1,1)=-33; % trial start well
            crun(1,3)=epochstart; % trial start time
            crun(1,2)=Results(firstreward,2); % trial end well
            crun(1,4)=Results(firstreward,1); % trial end time
            crun(1,5)=crun(1,4)-crun(1,3);
            
            for ii=size(crun,1)
                
                TS=crun(ii,3);
                TE=crun(ii,4);
                
                % get barrier on times
                if ~isempty(barrierevents);
                    BS=barrierevents(:,1);
                    BE=barrierevents(:,1)+barrierevents(:,4);
                    barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
                    crun(ii,6)=sum(barriertest)>0;
                else crun(ii,6)=0;
                end
            end
            
            run=crun;
            
            
        end
        % Find first trigger
        
        k=find(Results(:,3)>=0); % index of reward well triggers
        
        i=firstreward; % first rewarded index
        
        % Go through every index in k and get the corresponding times from
        % Results
        
        for i=i:(size(k,1)-1);
            crun=[];
            crun(1,1)=Results(k(i),2); % trial start well
            crun(1,3)=Results(k(i),5); % trial start time
            crun(1,2)=Results(k(i+1),2); % trial end well
            crun(1,4)=Results(k(i+1),1); % trial end time
            crun(1,5)=crun(1,4)-crun(1,3);
            
            % test if run contains a barrier
            % looks for overlap between run time and barrier times
            % get run start and run end times
            TS=crun(1,3);
            TE=crun(1,4);
            
            % get barrier on times
            if ~isempty(barrierevents);
                BS=barrierevents(:,1);
                BE=barrierevents(:,1)+barrierevents(:,4);
                barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
                crun(1,6)=sum(barriertest)>0;
            else crun(1,6)=0;
            end
            run=[run;crun];
            
            i=i+1;
            
        end
        
        if isempty(i)
            i=1;
        end
        
        % last trial
        drun=zeros(1,6);
        drun(1,1)=run(end,2); % trial start well
        drun(1,3)=Results(k(i),5); % trial start time
        drun(1,2)=-33; % trial end well
        drun(1,4)=epochend; % trial end time
        drun(1,5)=drun(1,4)-drun(1,3);
        % test if run contains a barrier
        % looks for overlap between run time and barrier times
        % get run start and run end times
        TS=drun(1,3);
        TE=drun(1,4);
        
        % get barrier on times
        if ~isempty(barrierevents);
            BS=barrierevents(:,1);
            BE=barrierevents(:,1)+barrierevents(:,4);
            barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
            drun(1,6)=sum(barriertest)>0;
        else drun(1,6)=0;
        end
        
        
        
        
        run=[run;drun];
        
        
        
    else
        % if no rewards are triggered, then the whole epoch is one run with
        % no defined wells
        crun=[];
        crun(1,1)=-33; % trial start well
        crun(1,3)=epochstart; % trial start time
        crun(1,2)=-33; % trial end well
        crun(1,4)=epochend; % trial end time
        crun(1,5)=crun(1,4)-crun(1,3);
        % test if run contains a barrier
        % looks for overlap between run time and barrier times
        % get run start and run end times
        TS=crun(1,3);
        TE=crun(1,4);
        
        % get barrier on times
        if ~isempty(barrierevents);
            BS=barrierevents(:,1);
            BE=barrierevents(:,1)+barrierevents(:,4);
            barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
            crun(1,6)=sum(barriertest)>0;
        else crun(1,6)=0;
        end
        
        
        
        
        
        run=crun;
        
    end
    
    
    
    
    
    % Make text version
    % Find start and end locations for each run
    if ~isempty(run)
        for i=3:5;
            Runtxt(:,i)=timetrans(run(:,i),10000,1);
        end
        
        for i=1:2;
            for j=1:size(run,1);
                [a,b]=find(activewells(:,2)==run(j,i));
                if size(a,1)~=0;
                    Runtxt{j,i}=strcat('Well',num2str(activewells(a,1)));
                else Runtxt{j,i}='-';
                end
                j=j+1;
            end
            i=i+i;
        end
        
        
        for j=1:size(run,1);
            if run(j,6)==1;
                Runtxt{j,6}='Barrier';
            else Runtxt{j,6}='-';
            end
            j=j+1;
        end
    else Runtxt=[];
    end
    
    
    % Write out stats
    
    Epochstats=[];
    epoch=epoch;
    duration=Results(size(Results,1),1)-Results(1,1);
    
    
    
    activewellinfo=[];
    
    for i=1:size(activewells,1);
        activewelltriggern=size(find(Wellevents(:,2)==activewells(i,2)),1);
        currwell(1,1)=activewells(i,1);
        currwell(1,2)=activewelltriggern;
        activewellrewardsn=size(find(Wellevents(:,3)==activewells(i,3)),1);
        currwell(1,3)=activewellrewardsn;
        activewellinfo=[activewellinfo;currwell];
        i=i+1;
    end
    
    rewardedwells=rewardedwells;
    
    
else % for sleep sessions
    
    duration=ranges(epoch,2)-ranges(epoch,1);
    Config=[];
    Configinfo=[];
    Wellevents=[];
    barrierevents=[];
    outputtriggertimes=[];
    
    Results=[];
    ResultsTxt=[];
    Runtxt=[];
    activewells=[];
    rewardedwells=[];
    run=[];
    epochtype='Sleep';
    
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% Use Posepochprocess to write out position data

daydirect=strcat(datadir,animalname,'/',animalname,dsz,dayt,'/');

day=strcat(dsz,dayt);

posout=Posepochprocess_realigntimestamp(daydirect,animalname,day,epoch,'cmperpix',pixdim);

% output from Posepochprocess should contain both the rawpos and the
% corrected pos data

% lens correction done in posinterpcorrection
% % lens distorsion correction
%
% rawpos(:,1)=pos.data(:,2)/pixdim;
% rawpos(:,2)=pos.data(:,3)/pixdim;
%
% % converts into indices for distorsion correction
%
% pointind=pixelindex(rawpos,320,240);
%
% % lookup coordinates of the transformed indices
%
% coord=lenscorrtransf(pointind);
%
% % put results in new corrected data matrix
%
% pos.correcteddata(:,1)=pos.data(:,1);
% pos.correcteddata(:,2)=coord(:,1)*0.5;
% pos.correcteddata(:,3)=coord(:,2)*0.5;
% pos.correcteddata(:,4)=pos.data(:,4);
% pos.correcteddata(:,5)=pos.data(:,5);

% position data after lens distortion correction
pos.correcteddata=posout.pos.data;
pos.correcteddatainfo=['Lens distortion correction in home/Src/ImageCorrection'];
pos.arg=posout.pos.arg;
pos.descript=posout.pos.descript;
pos.fields=posout.pos.fields;
pos.cmperpixel=posout.pos.cmperpixel;
pos.rawpos=posout.rawpos.data; % raw position data in pixel dimensions
pos.rawposdescript=['timestamp x1 y1 x2 y2 - in pixels'];
pos.headdirection=posout.pos.direction;

%pos.SKpos=posout.SKpos;

Events={};

%Epochstats=[animal;epoch;duration];
Epochdata.Config=Config;
Epochdata.Configinfo=Configinfo;
Epochdata.Events.Wellevents=Wellevents;
Epochdata.Events.Barrier=barrierevents;
Epochdata.Events.Welltriggers=outputtriggertimes;
Epochdata.Stats.animal=epochanimal;
Epochdata.Stats.epoch=epoch;
Epochdata.Stats.epochtype=epochtype;
Epochdata.Stats.duration=duration;

Epochdata.Results=Results;
Epochdata.Resultstxt=ResultsTxt;
Epochdata.Runtxt=Runtxt;
Epochdata.Runtxtinfo={'Start','Goal','Start time','End time','Duration','Barrier'};
Epochdata.Wellinfo.rewardedwells=rewardedwells;

Epochdata.Run=run;
Epochdata.Pos=pos;

fileoutdir=strcat(datadir,animaldir,'/');

end


% change to that directory and saves the figure with file name
% animal_day_epoch

%cd(fileoutdir);
%filename = strcat(animalname,'_',dayt,'e',epocht,'.mat');
%save(filename,'Epochdata');


%cd(currentdir);






