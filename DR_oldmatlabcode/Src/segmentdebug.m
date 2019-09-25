%% find segments used during each trial
%code taken from PlotRerouteSeg_multiday
%used to break down each trial into individual segments
%Trialsegments is the sequence of segments used in each trial
%Trialsegments{trial}([time startintersection endintersection segment])

%open segmenttable, identifies a track segment from two
%track intersections

load('/home/jai/Src/NSpikeProcess/segmentlookuptable.mat');

day=9;
epoch=6;
index(1)=day;
index(2)=epoch;

dsz = '';
if (index(1) < 10)
    dsz = '0';
end

directoryname='/data14/jai/M2_/';
fileprefix='M2';



%% load the data and extract position values
load(strcat(directoryname,fileprefix,'data', dsz, num2str(index(1)), '.mat'));
%load(strcat(directoryname,fileprefix,'pos03', '.mat'));
%eval(['Data = ',lowercasethree,'pos;'])
load(strcat(directoryname,fileprefix,'task',dsz, num2str(index(1)), '.mat'));
load(strcat(directoryname,fileprefix,'linpos',dsz, num2str(index(1)), '.mat'));

trajmatrix=linpos{index(1)}{index(2)}.trajmatrix;

% re done to extract start and end of each trial from data

%rewardedevents=Data{1,day}{1,epoch}.Events.Welltriggers; % reward events as reward release rather than reward trigger
rewardedevents=data{1,day}{1,epoch}.Run;

if ~isempty(rewardedevents);
    %rewardedevents=rewardedevents((rewardedevents(:,4)>100),:);
    Trialsegments=[]; % table containing the start time and end time of a trial and sequence of segments used in that trial
    barrierruns=[];
    % get start and end times for each run
    %             for i=1:size(Data{1,daylist(1,d)}{1,epoch}.Run,1);
    %                 tstart=Data{1,daylist(1,d)}{1,epoch}.Run(i,3);
    %                 tend=Data{1,daylist(1,d)}{1,epoch}.Run(i,4);
    %
    for i=1:size(rewardedevents,1);
        
        timeind=rewardedevents(i,3);
        tstartref=min(find(data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000>(timeind))); % for some reason +10000 seems to be accurate rather than no offset
        %timeind2=rewardedevents(i+1,1)+10000; % when rat first triggers reward well, location is not within defined well region, so take location at 10000 later
        timeind2=rewardedevents(i,4);
        tendref=max(find(data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000<(timeind2)));
        % since all segment indices are not in time units,
        % need to convert time back to indices.
        
        tstartind=data{1,day}{1,epoch}.Pos.correcteddata(tstartref,1)*10000;
        tendind=data{1,day}{1,epoch}.Pos.correcteddata(tendref,1)*10000;
        
        %seg=segment(tstartref:tendref,1);
        % get exit times during period
        
        estartref=find(trajmatrix(:,1)*10000>=timeind); %removed +10000 since the epoch times are more reliablely calculated in batchDIOrerouteextract
        if ~isempty(estartref);
            estartref=estartref(1,1);
            [val eendref]=min(abs(trajmatrix(:,1)*10000-timeind2));
            if ~isempty(estartref);
                exitref=trajmatrix(estartref(1,1):eendref,:);
            else i=size(data{1,day}{1,epoch}.Run,1);
            end
            % concatenate columns to generate integer representation of
            % path, this is the sequence of segments used
            % in a trial
            for q=1:size(exitref,1);
                exitref(q,4)=str2num(sprintf('%d%d',exitref(q,2),exitref(q,3)));
            end
            %                             % remove segments containing
            %                             barrier
            %                             for b=1:size(barrier,1);
            %                                 if barrier(b,1)==daylist(1,d);
            %                                     if barrier(b,2)==epoch;
            %                                         fdirect=str2num(sprintf('%d%d',barrier(b,3),barrier(b,3)));
            %                                         forward=find(exitref(:,4)~=fdirect);
            %                                         exitref=exitref(forward,:);
            %                                         rdirect=str2num(sprintf('%d%d',barrier(b,4),barrier(b,4)));
            %                                         reverse=find(exitref(:,4)~=rdirect);
            %                                         exitref=exitref(reverse,:);
            %                                     end
            %                                 end
            %                                 b=b+1;
            %                             end
            %
            % remove incomplete path crossings
            exitref=exitref(find(exitref(:,2)-exitref(:,3)~=0),:);
            
            % find segment number from connectivity matrix rcon
            for q=1:size(exitref,1);
                exitref(q,5)=segmentlookuptable(rowfind(exitref(q,4),segmentlookuptable(:,1)),2);
                %exitref(q,6)=disttb(exitref(q,5),1);
            end
            % remove incomplete path crossings
            exitref=exitref(find(exitref(:,2)-exitref(:,3)~=0),:);
            
            % removes repeated use of corner segments
%             exclude_seg = [9 10 11 12];
%             exitref_tmp = exitref(:,5);
%             exitref_tmp(~ismember(exitref(:,5),exclude_seg)) = repmat(0,1,length(find(~ismember(exitref(:,5),exclude_seg))));
%             if(length(exitref_tmp) > 1)
%                 exitref(find(diff(exitref_tmp)==0 & exitref_tmp(1:end-1) ~= 0),:) = [];
%             end
            
            Trialsegments{1,i}=exitref(:,[1 2 3 5]); % put values in matrix
            %segments(i,2)=sum(exitref(:,6)); % distance for path
            
            % get barrier index for each run
            % calculate run by run if barrier is on
            
            %                             barrierevents=Data{1,day}{1,epoch}.Events.Barrier;
            %                             if ~isempty(barrierevents);
            %
            %                                 % get run start and run end times
            %                                 TS=Data{1,day}{1,epoch}.Pos.correcteddata(tstartref,1)*10000;
            %                                 TE=Data{1,day}{1,epoch}.Pos.correcteddata(tendref,1)*10000;
            %                                 BS=barrierevents(:,1);
            %                                 BE=barrierevents(:,1)+barrierevents(:,4);
            %                                 barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
            %                                 if sum(barriertest)>0;
            %                                     segments(i,3)=1;
            %                                 else segments(i,3)=0;
            %                                 end
            %                             else segments(i,3)=0;
            %                             end
            
            
            if isempty(Trialsegments{1,i});
                display(sprintf('Respecify well location for day %d epoch %d reward index %d start at %d end at %d',day,epoch,i,tstartref,tendref));
            end
            i=i+1;
        end
        
    end
    %else k=k+1;
end