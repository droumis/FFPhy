function [] = PlotRerouteSeg_multiday( animaldir,animalname,daylist, barrier)
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
% needs segment table located in
%/home/jai/Src/NSpikeProcess/segmenttable.mat

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
% stats for barriers for each run of each day
daybarrierruns=[];
% load connectivity matrix and distance table
load('/home/jai/Src/NSpikeProcess/rerouteconnectivity.mat');
rcon=rerouteconnectivity;
clear rerouteconnectivity;
% load segmenttable matrix and distance table
load('/home/jai/Src/NSpikeProcess/segmenttable.mat');
for d=1:size(daylist,2);
    % open files
    dayt = num2str(daylist(1,d));
    % ---- File loading ----
    % See if day number needs 0
    dsz = '';
    if (daylist(1,d) < 10)
        dsz = '0';
    end
    % Load the mat file
    datadir = '/data14/jai/'; % Specify data directory and load the file
    dayt = num2str(daylist(1,d)); % Converts the day and epoch into text
    sfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'_',dsz,dayt,'.mat');
    load(sfilename);
    posfilename = strcat(datadir,animaldir,'/',animaldir(1:end-1),'linpos',dsz,dayt,'.mat');
    load(posfilename);
    % find no. epochs
    j =  size(Data{1,daylist(1,d)},2);
    % data for all epochs of each day
    epochrunmat=[];
    epochmat=[];
    epochstat=[];
    currdaybarrierruns=[];
    segusage=[];
    for k=1:size(Data{1,daylist(1,d)},2);
        if isempty(Data{1,daylist(1,d)}{1,k});
            k=k+1;
        elseif strcmp(Data{1,daylist(1,d)}{1,k}.Stats.epochtype,'Sleep')==1;
            k=k+1;
        else
            if strcmp(Data{1,daylist(1,d)}{1,k}.Stats.animal(1,1),animalname)==1; %find all epochs for animal
                epoch = k;
                epocht = num2str(epoch);
                %rewardedevents=Data{1,day}{1,epoch}.Events.Wellevents((Data{1,day}{1,epoch}.Events.Wellevents(:,3)~=-50),:);
                
                rewardedevents=Data{1,day}{1,epoch}.Events.Welltriggers; % reward events as reward release rather than reward trigger
                if ~isempty(rewardedevents);
                    rewardedevents=rewardedevents((rewardedevents(:,4)>100),:);
                    segments=[]; %clean matrix
                    barrierruns=[];
                    % get start and end times for each run
                    %             for i=1:size(Data{1,daylist(1,d)}{1,epoch}.Run,1);
                    %                 tstart=Data{1,daylist(1,d)}{1,epoch}.Run(i,3);
                    %                 tend=Data{1,daylist(1,d)}{1,epoch}.Run(i,4);
                    %
                    for i=1:size(rewardedevents,1)-1;
                        timeind=rewardedevents(i,1)+rewardedevents(i,4);
                        tstartref=min(find(Data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000>(timeind))); % for some reason +10000 seems to be accurate rather than no offset
                        timeind2=rewardedevents(i+1,1)+10000; % when rat first triggers reward well, location is not within defined well region, so take location at 10000 later
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
                            segusage=[segusage;exitref(:,4)];
                            
                            
                            % get barrier index for each run
                            % calculate run by run if barrier is on
                            
                            barrierevents=Data{1,day}{1,epoch}.Events.Barrier;
                            if ~isempty(barrierevents);
                                
                                % get run start and run end times
                                TS=Data{1,day}{1,epoch}.Pos.correcteddata(tstartref,1)*10000;
                                TE=Data{1,day}{1,epoch}.Pos.correcteddata(tendref,1)*10000;
                                BS=barrierevents(:,1);
                                BE=barrierevents(:,1)+barrierevents(:,4);
                                barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
                                if sum(barriertest)>0;
                                    segments(i,3)=1;
                                else segments(i,3)=0;
                                end
                            else segments(i,3)=0;
                            end
                            
                            
                            if segments(i,1)==0;
                                display(sprintf('Respecify well location for day %d epoch %d reward index %d start at %d end at %d',day,epoch,i,tstartref,tendref));
                            end
                            i=i+1;
                        end
                    end
                    
                    %currbarrierruns=Data{1,daylist(1,d)}{1,epoch}.Run(:,6);
                    %barrierruns=[barrierruns;currbarrierruns];
                    epochrunmat=[epochrunmat;size(Data{1,daylist(1,d)}{1,epoch}.Run,1)];
                    epochmat=[epochmat;segments];
                    
                    % matrix with [day epoch run_number]
                    currepochstat=[];
                    currepochstat(1,1)=daylist(1,d);
                    currepochstat(1,2)=epoch;
                    currepochstat(1,3)=size(Data{1,daylist(1,d)}{1,epoch}.Run,1);
                    epochstat=[epochstat;currepochstat];
                    %currdaybarrierruns=[currdaybarrierruns;barrierruns];
                    k=k+1;
                
                end
            else display(sprintf('No %s data from day %s epoch %s',animalname,num2str(daylist(1,d)),num2str(k)));
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
        daybarrierruns=[daybarrierruns;currdaybarrierruns];
        d=d+1;
    end


% calculates number of times each segment is used, directional

[usage, segmentind]=hist(segusage,sort(segmenttable(:,1)));
%[usage, segmentind]=hist(segusage,unique(segusage));
segmentusage=[segmentind usage'];
% segmentusage=zeros(max(segmentind),3);
% segmentusage(segmentind)=segmentusagetemp((segmentusagetemp(:,1)==segmentind),2);

for i=1:size(segmentusage);
    segmentusage(i,3)=segmenttable(segmenttable(:,1)==segmentusage(i,1),2);
    segmentusage(i,4)=segmenttable(segmenttable(:,1)==segmentusage(i,1),3);
end

% segmentusage(segmentind,2)=segmenttable((segmenttable(:,1)==segmentind),2);
% segmentusage(segmentind,2)=segmenttable((segmenttable(:,1)==segmentind),3);
% segmenttable(:,4)=usage(segment(:,1)==segmenttable(:,1));
segmentusage=sortrows(segmentusage,[3,4]);

figure;

bar([segmentusage(segmentusage(:,4)==1,2) segmentusage(segmentusage(:,4)==0,2)],'stacked');
xlabel('Path segment');
ylabel('Number of passes');
legend('Forward','Reverse', 'Location','EastOutside');
title(sprintf('Path usage of %s on day %s',animalname,dayt));
d=d+1;
end
end


% ------------------



