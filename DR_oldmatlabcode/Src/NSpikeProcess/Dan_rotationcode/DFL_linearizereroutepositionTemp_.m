% [statematrix, segmenttable, trajwells,wellSegmentInfo, segmentInfo,trajmatrix] = LINEARIZEPOSITION(directoryname,fileprefix,index, options)
%
% Takes each position from the animal's 'pos' structure  and
% projects it onto the linear trajectory segments for the environment
% and task defined by coordprogram.
% This routine finds the segment closest to each position, projects
% the position to the closest location on the segment. It then calculates
% the linear distance from each well to the animal and which well-to-well
% trajectory the animal is on for each time step.
%
%      --INPUTS--
%
%    directoryname - example 'data99/user/animaldatafolder/', a folder
%                    containing processed matlab data for the animal
%       fileprefix - animal specific prefix for each datafile
%            index - [day epoch]
%
%          options - 'lowercasethree' - varaible prefix (default '')
%
%                    'maxvelocity' - maximum allowed velocity in cm/sec (default 300)
%
%                    'maxsegdiff' - maximum allowed number of segments that
%                    could be traversed from one time point to the next.  This
%                    is useful primarily for linear environments where there is
%                    no branching. Default 100
%
%                    'welldist' - to decide whether or not the animal
%                    completed the trajectory, the program needs a radius
%                    around the endpoint of each trajectory which the
%                    animal must enter to trigger a completed trajectory (default 5 cm)
%
%                    'mindiff' - the minimum amount of time between two
%                    foodwell zone entries to be concidered a real
%                    trajectory (default 2 seconds)
%
%                    'velocitysmoothwindow' - to compute linear speed, the
%                    program needs to smooth linear position data.  It is smoothed
%                    with a gaussian of length VSW and std VSW/4 (default 2
%                    seconds)
%
%                    'welldistthresh' -the threshold linear distance (in cm) from the well for the program
%                    to add a trajectory when the animal turns around
%                    before reaching the end well.  (default 0)
%
%                   'branchpointdist' -this performs a similar function to
%                   welldistthresh, except the distance is defined as the
%                   distance traveled a the part of the track that only leads to one well.
%                   If a value is given for this, it supercedes any value
%                   given for welldistthresh.
%
%
%       --OUTPUTS--
%
%      statematrix - a structure with each field the same length as pos.data
%                    Contains information about the linear position of the anuimal
%                    for each time step.
%
%     segmenttable - a lookup table with 3 columns: the first column
%                    is the linear trajectory number (the row number in
%                    lindistpos{i}{e}.traject), the second column is the segment number in
%                    that trajectory, and the third column is the unique indentifier number
%                    for that segment, which is used in the segmentIndex field of statematrix. Any segments
%                    that are used in multiple trajectories are given the same identifier number.
%
%        trajwells - gives the start and end well numbers for each trajectory (each row is
%                    a trajectory).  Each trajecory is defined from the
%                    output of coordprogram.
%

function [statematrix, segmenttable, trajwells, wellSegmentInfo, segmentInfo, ...
    trajmatrix, Trialsegments, Trialsegments_raw] = DFL_linearizereroutepositionTemp(directoryname,fileprefix, index, varargin)

% if rewarded events is empty, make sure all return values still return
Trialsegments = [];
Trialsegments_raw = []; % for debugging segment removal, see later

maxv = 1000;
maxsegdiff = 100;
lowercasethree = ''; %default variable prefix is none
welldist = 10;
mindiff = 0.5;
welldistthresh = 0;
branchpointdist = inf;
smoothwidth = 1;
if (length(varargin) == 1)
    varargin = varargin{1};
end
%set variable options
for option = 1:2:length(varargin)-1
    
    switch varargin{option}
        case 'lowercasethree'
            lowercasethree = varargin{option+1};
        case 'maxvelocity'
            maxv = varargin{option+1};
        case 'maxsegdiff'
            maxsegdiff = varargin{option+1};
        case 'welldist'
            welldist = varargin{option+1};
        case 'mindiff'
            mindiff = varargin{option+1};
        case 'velocitysmoothwindow'
            smoothwidth = varargin{option+1};
        case 'welldistthresh'
            welldistthresh = varargin{option+1};
        case 'branchpointdist'
            branchpointdist = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

dsz = '';
if (index(1) < 10)
    dsz = '0';
end
%% load the data and extract position values
load(strcat(directoryname,fileprefix,'data', dsz, num2str(index(1)), '.mat'));
%load(strcat(directoryname,fileprefix,'pos03', '.mat'));
%eval(['Data = ',lowercasethree,'pos;'])
load(strcat(directoryname,fileprefix,'task',dsz, num2str(index(1)), '.mat'));
%eval(['task = ',lowercasethree,'task;'])
Pos = data{index(1)}{index(2)}.Pos;
pos = Pos.correcteddata;
task = task{index(1)}{index(2)};
timestep = pos(2,1) - pos(1,1);
% get reward wells
wellno=data{1,index(1)}{1,index(2)}.Wellinfo.rewardedwells;

%% get track type, specific to reroute
track=task.type;
% make sure we have the direction information
toknum = isdatafield(Pos.fields, 'dir');
%if (~toknum)
%    Error('No direction field in pos');
%end

%% initialize variables
poslength = size(pos,1);
segment = zeros(poslength,1);
lindist = ones(size(pos,1), 1) * -1;
vect = zeros(poslength,2);
newpos = pos;

%% use the coordinate program to fetch the track coordinates
% array has task{day}{epoch}.linearcoord{segment no.}{coordinates of segment ends at all pos}
coordInTime = task.linearcoord;
clear task;

%eval(['coordProgramHandle = @',coordprogram,';']);
%coordInTime = feval(coordProgramHandle,pos(:,1),task);

%% calculate the linear distance for each segment [0;distance]
distsum = [];
for i = 1:length(coordInTime)
    firstcoord{i} = coordInTime{i}(:,:,1);
    if (size(firstcoord{i},1) > 1)
        distsum{i}(1,1) = 0;
        for j = 2:size(firstcoord{i},1)
            distsum{i}(j,1) = sqrt( ((firstcoord{i}(j,1) - firstcoord{i}(j-1,1))^2) + ((firstcoord{i}(j,2) - firstcoord{i}(j-1,2))^2) ) + ...
                distsum{i}(j-1);
        end
    end
end

%% Process trajectory data
%create a table of segments describing which trajectory and segnum
% they are.  This table give segment an identifying number. If any segment
%is used in multiple trajectories, it is only listed once as belonging to
%the first trajectory. This function also calculates specific information
%about the track segments and trajectories.
[segmenttable, intersectcoord, wellindex, trajwells, wellSegmentInfo, segmentInfo] = getSegmentTable(firstcoord,track,wellno);

ntraj = size(trajwells,1);

%get the well locations for all time frames
for i = 1:size(wellindex,1)
    welllocations(i,1:2,:) = coordInTime{wellindex(i,1)}(wellindex(i,2),1:2,:);
end



% the main loop projects each position point onto each segment (and does error correction along the way)
inbound = 0;
lastvalid = -1;
lastsegind = -1;
for i = 1:size(pos,1)
    %lookup the segment coordinates for this time point
    for findcoord = 1:length(segmentInfo.segmentLength)
        tableInd = min(find(segmenttable(:,3) == findcoord)); % index of the segment
        coord{findcoord} = coordInTime{segmenttable(tableInd,1)}(segmenttable(tableInd,2):segmenttable(tableInd,2)+1,:,i); % coordinates of the end points
        %calculate each segment's vector
        coordvector(findcoord,1:2) = diff(coord{findcoord});
    end
    
    if ~(mod(i,1000))
        %display percent progress
        disp([num2str(round((i/poslength)*100)),'%']);
    end
    % project each position point onto the linear segments
    if (pos(i,2) == 0)
        % this is an invalid point, so set newpos to zero for this position
        newpos(i, 2:3) = 0;
        warning(sprintf('zero position at index %d, animal %s, day epoch %d %d', i, fileprefix, index(1), index(2)));
    else
        tmppos = [];
        for j = 1:length(coord)
            %tmppos(j,1:5) = projectpoint(pos(i,2:3), coord{j}); %a 1x5 vector with x, y, distance, onseg, and segnum
       tmppos(j,1:5) = JY_projectpoint(pos(i,2:3), coord{j}); %a 1x5 vector with x, y, distance, onseg, and segnum
        
        end
        %change the segment number to equal the segment index
        tmppos(:,5) = 1:size(tmppos,1);
        % take the point that is the least distance from the segments and that
        % does not force an excessive velocity from the last point
        
        tmppos = sortrows(tmppos,3);
        % check the velocity to the last valid point
        
        if (lastvalid ~= -1)
            tmppnt = [0 0];
            tmppnt(1,1) = newpos(lastvalid,2);
            tmppnt(1,2) = newpos(lastvalid,3);
            vel = dist(tmppnt,tmppos(:,1:2)) ./ (newpos(i,1) - ...
                newpos(lastvalid,1));
            % get the number of segments we will have moved across
            segdiff = abs(tmppos(:,5) - lastsegind);
            % the next point is the first point in the list with a
            % corresponding velocity less than maxv and a change in the
            % number of segments
            % less than maxsegdiff
            newind = min(find((vel < maxv) & (segdiff <= maxsegdiff)));
        else
            %no last valid segment has yet been found
            newind = 1;
        end
        
        if (isempty(newind))
            % set this point to zero, as it is not a valid point
            disp(['Warning: undetermined linear position at pos index ',num2str(i)]);
            newpos(i,2:3) = [0 0];
        else
            newsegment = tmppos(newind,5);
            segment(i) = newsegment;
            segdist(i) = sqrt(sum((tmppos(newind,1:2) - coord{newsegment}(1,:)).^2)); %the distance along the segment
            vect(i,1) = coordvector(newsegment,1);
            vect(i,2) = coordvector(newsegment,2);
            newpos(i,2:3) = round(tmppos(newind, 1:2));
            
            %find the linear distance for the point from each well
            %this requires a check for which direction the segment is
            %aligned relative to the well
            
            %--------------------------------------------------------------
            
            
            if track~='reroute';
                distToSegment = wellSegmentInfo.distanceTable(:,newsegment)';
                segmentDist = [];
                for wellcount = 1:size(wellSegmentInfo.distanceTable,1)
                    if  (wellSegmentInfo.segmentDirection(wellcount,newsegment) == 1)
                        %it is aligned in the foreward direction
                        segmentdist(1,wellcount)  = segdist(i);
                    else
                        %it is aligned in the backward direction
                        segmentdist(1,wellcount)  = segmentInfo.segmentLength(newsegment)-segdist(i);
                    end
                end
                
                %calculate the linear distance from each well
                lindist(i,1:wellcount) = distToSegment + segmentdist;
            end
            %--------------------------------------------------------------
            lastvalid = i;
            lastsegind = segment(i);
        end
        
    end
    
    
end

%  fill in any missing elements skipped from invalid positions

segdist = vectorfill(segdist,-1);
vect(:,1) = vectorfill(vect(:,1),0);
vect(:,2) = vectorfill(vect(:,2),0);
segment = vectorfill(segment,0);
for i = 1:size(lindist,2)
    lindist(:,i) = vectorfill(lindist(:,i),-1);
end

%% Getting trajectories of individual trials
% segment at each point is specified by 'segment'
% use this information to find sequence of segments in each trial



%% get the foodwell enter and exit times: [starttime startwell endwell]
[trajmatrix trajmatrix_raw] = gettraject(pos, welllocations,trajwells, lindist, wellSegmentInfo, welldist, mindiff, welldistthresh,branchpointdist,track,intersectcoord);



%% find segments used during each trial
%code taken from PlotRerouteSeg_multiday
%used to break down each trial into individual segments
%Trialsegments is the sequence of segments used in each trial
%Trialsegments{trial}([time startintersection endintersection segment])

%open segmenttable, identifies a track segment from two
%track intersections

load('/home/jai/Src/NSpikeProcess/segmentlookuptable.mat');



day=index(1);
epoch=index(2);

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
        
        seg=segment(tstartref:tendref,1);
        % get exit times during period
        
        estartref=find(trajmatrix(:,1)*10000>=timeind); %removed +10000 since the epoch times are more reliablely calculated in batchDIOrerouteextract
        if ~isempty(estartref);
            estartref=estartref(1,1);
            diffvals=trajmatrix(:,1)*10000-timeind2;
            eendref=max(find(diffvals<=0));
            if isempty(eendref)
                [vals eendref]=max(diffvals);
            end
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
%             exitref_tmp(~intersect(exitref(:,5),exclude_seg)) = repmat(0,1,length(find(~intersect(exitref(:,5),exclude_seg))));
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

% this section is for debugging segment assignment
% trialsegments_raw is with repeated segment usage kept

% rewardedevents=Data{1,day}{1,epoch}.Events.Welltriggers; % reward events as reward release rather than reward trigger
% if ~isempty(rewardedevents);
%     rewardedevents=rewardedevents((rewardedevents(:,4)>100),:);
%     Trialsegments_raw=[]; % table containing the start time and end time of a trial and sequence of segments used in that trial
%     barrierruns_raw=[];
%     % get start and end times for each run
%     %             for i=1:size(Data{1,daylist(1,d)}{1,epoch}.Run,1);
%     %                 tstart=Data{1,daylist(1,d)}{1,epoch}.Run(i,3);
%     %                 tend=Data{1,daylist(1,d)}{1,epoch}.Run(i,4);
%     %
%     for i=1:size(rewardedevents,1)-1;
%         timeind=rewardedevents(i,1)+rewardedevents(i,4);
%         tstartref=min(find(Data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000>(timeind))); % for some reason +10000 seems to be accurate rather than no offset
%         %timeind2=rewardedevents(i+1,1)+10000; % when rat first triggers reward well, location is not within defined well region, so take location at 10000 later
%         timeind2=rewardedevents(i+1,1);
%         tendref=max(find(Data{1,day}{1,epoch}.Pos.correcteddata(:,1)*10000<(timeind2)));
%         % since all segment indices are not in time units,
%         % need to convert time back to indices.
%         
%         tstartind=Data{1,day}{1,epoch}.Pos.correcteddata(tstartref,1)*10000;
%         tendind=Data{1,day}{1,epoch}.Pos.correcteddata(tendref,1)*10000;
%         
%         seg=segment(tstartref:tendref,1);
%         % get exit times during period
%         
%         estartref=find(trajmatrix(:,1)*10000>timeind+10000); %removed +10000
%         if ~isempty(estartref);
%             estartref=estartref(1,1);
%             [val eendref]=min(abs(trajmatrix(:,1)*10000-timeind2));
%             if ~isempty(estartref);
%                 exitref=trajmatrix(estartref(1,1):eendref,:);
%             else i=size(Data{1,day}{1,epoch}.Run,1);
%             end
%             % concatenate columns to generate integer representation of
%             % path, this is the sequence of segments used
%             % in a trial
%             for q=1:size(exitref,1);
%                 exitref(q,4)=str2num(sprintf('%d%d',exitref(q,2),exitref(q,3)));
%             end
%             %                             % remove segments containing
%             %                             barrier
%             %                             for b=1:size(barrier,1);
%             %                                 if barrier(b,1)==daylist(1,d);
%             %                                     if barrier(b,2)==epoch;
%             %                                         fdirect=str2num(sprintf('%d%d',barrier(b,3),barrier(b,3)));
%             %                                         forward=find(exitref(:,4)~=fdirect);
%             %                                         exitref=exitref(forward,:);
%             %                                         rdirect=str2num(sprintf('%d%d',barrier(b,4),barrier(b,4)));
%             %                                         reverse=find(exitref(:,4)~=rdirect);
%             %                                         exitref=exitref(reverse,:);
%             %                                     end
%             %                                 end
%             %                                 b=b+1;
%             %                             end
%             %
%             % remove incomplete path crossings
%             exitref=exitref(find(exitref(:,2)-exitref(:,3)~=0),:);
%             
%             % find segment number from connectivity matrix rcon
%             for q=1:size(exitref,1);
%                 exitref(q,5)=segmentlookuptable(rowfind(exitref(q,4),segmentlookuptable(:,1)),2);
%                 %exitref(q,6)=disttb(exitref(q,5),1);
%             end
%             % remove incomplete path crossings
%             exitref=exitref(find(exitref(:,2)-exitref(:,3)~=0),:);
%             
%             exclude_seg = [9 10 11 12];
%             exitref_tmp = intersect(exitref(:,5),exclude_seg) & exitref(:,5);
%             exitref(find(diff(exitref_tmp)==0),:) = [];
%             
%             Trialsegments_raw{1,i}=exitref(:,[1 2 3 5]); % put values in matrix
%             %segments(i,2)=sum(exitref(:,6)); % distance for path
%             
%             % get barrier index for each run
%             % calculate run by run if barrier is on
%             
%             %                             barrierevents=Data{1,day}{1,epoch}.Events.Barrier;
%             %                             if ~isempty(barrierevents);
%             %
%             %                                 % get run start and run end times
%             %                                 TS=Data{1,day}{1,epoch}.Pos.correcteddata(tstartref,1)*10000;
%             %                                 TE=Data{1,day}{1,epoch}.Pos.correcteddata(tendref,1)*10000;
%             %                                 BS=barrierevents(:,1);
%             %                                 BE=barrierevents(:,1)+barrierevents(:,4);
%             %                                 barriertest=(TS >= BS & TS < BE) | (TE <= BE & TE > BS) | (TS <= BS & TE >= BE);
%             %                                 if sum(barriertest)>0;
%             %                                     segments(i,3)=1;
%             %                                 else segments(i,3)=0;
%             %                                 end
%             %                             else segments(i,3)=0;
%             %                             end
%             
%             
%             if isempty(Trialsegments_raw{1,i});
%                 display(sprintf('Respecify well location for day %d epoch %d reward index %d start at %d end at %d',day,epoch,i,tstartref,tendref));
%             end
%             i=i+1;
%         end
%         
%     end
%     %else k=k+1;
% end

%Now, we want to find out which well to use as the reference point for
%linear distance traveled.  This depends on which trajectory the animal is
%on.
whichTraj = zeros(size(pos,1),1);

%----------------------------------------------------------------------
%if track~='reroute';

%if all trajectories start with the same well, then we should always use
%that well as the reference point.
if ( (sum(abs(diff(trajwells(:,1)))) == 0)| (isempty(diff(trajwells(:,1)))) )
    whichTraj(:) = trajwells(1,1); %just pick the first traj for all points, later the reference well will be chosen as the trajectory's starting point
end
if ~isempty(trajmatrix)
    %for each timestamp in pos, get which of these complete trajectories the
    %animal is in
    welltraj = trajmatrix(lookup(pos(:,1),trajmatrix(:,1),-1),2:3);
else
    disp('WARNING: No completed trajectories.  Check endpoint locations if this is incorrect.')
    %because no trajectories were competed, we fill the start and end wells
    %with NaN's
    welltraj = nan(size(pos,1),2);
end
samewelltraj = [];

%end
%-----------------------------------------------------------------------

%
% if (whichTraj(1) == 0) %more than one starting well exists, so we need more complicated work to find the reference well
%
%     %if the animal was on a segment that belongs to only one trajectory,
%     %then we assign that trajectory to that time
%     singleTrajSegments = find(rowcount(segment,segmenttable(:,3))==1);
%     whichTraj(singleTrajSegments) = segmenttable(rowfind(segment(singleTrajSegments),segmenttable(:,3)),1);
%     zeroind = find(whichTraj == 0); %these are the indeces to the times when the animal was on amiguous segments
%
%     %next, we look at the times when the animal was on segments belonging
%     %to 2 or more trajectories.  We use the start and end wells for each
%     %time step to determine the best trajectory to assign the timestep to.
%
%     %-------------------------------------------------------------------
%     if track~='reroute';
%
%         if ~isempty(trajmatrix)
%
%             %find the trajectory number for each well exit-to-entry, if it exists.
%             %to to this we must check both the foreward and reverse directions
%             tmptrajnum = rowfind(welltraj(zeroind,:),trajwells);
%             tmptrajnum2 = rowfind(welltraj(zeroind,:),trajwells(:,[2 1]));
%             tmptrajnum = max(tmptrajnum,tmptrajnum2);
%             trajNonzeros = find(tmptrajnum > 0);
%             %only include the times when the segment belonged to the trajectory
%             segmenttraj = segmentInfo.segmentTraj(segment);
%             tmptrajnum(trajNonzeros) = tmptrajnum(trajNonzeros) .* bitget(segmenttraj(zeroind(trajNonzeros)),tmptrajnum(trajNonzeros));
%             whichTraj(zeroind) = tmptrajnum;
%             zeroind = find(whichTraj == 0);
%             samewelltraj = welltraj(zeroind(welltraj(zeroind,1) == welltraj(zeroind,2)),1);
%             samewelltrajind = zeroind(welltraj(zeroind,1) == welltraj(zeroind,2));
%             %only use the same-well trajectories where the well is found as the
%             %start point for one of the defined trajectories
%             foundstart = find(intersect(samewelltraj,trajwells(:,1)));
%             samewelltraj = samewelltraj(foundstart);
%             samewelltrajind = samewelltrajind(foundstart);
%         end
%     end
%     %------------------------------------------------------------------
% end

%% Creating referencewell list
%  referencewell is always set as the entering node of the current segment

referencewell = welltraj(:,1);

%% Calculate linpos relative to segment the rat is on with distance 0 fixed
%% on the first segment point listed in segmentInfo.segmentCoords
lindist = zeros(size(pos,1),1);
temp_count = 0;
for ii = 1:size(pos,1)
    %refwell_coord_temp = intersectcoord(referencewell(ii),:);
    refwell_coord_temp = segmentInfo.segmentCoords(segment(ii),1:2);
    pos_coord_temp = pos(ii,2:3);
    pos_vec_temp = pos_coord_temp - refwell_coord_temp;
    segment_coord_temp = segmentInfo.segmentCoords(segment(ii),:);
    segment_vec_temp = [segment_coord_temp(3) - segment_coord_temp(1) segment_coord_temp(4) - segment_coord_temp(2)];
    
    lindist(ii) = dot(pos_vec_temp,segment_vec_temp)/norm(segment_vec_temp);
    
    if(segmentInfo.segmentLength(segment(ii)) < lindist(ii) | segmentInfo.segmentLength(segment(ii)) < 0)
        temp_count = temp_count + 1;
        %disp(sprintf('pos %d (%f, %f), seg %d, seglength %f, lindist %f',ii,pos(ii,2),pos(ii,3),segment(ii),segmentInfo.segmentLength(segment(ii)),lindist(ii)));
        %disp(sprintf('pos_vec %s, segment_vec %s, refwell %s',mat2str(pos_vec_temp),mat2str(segment_vec_temp),mat2str(refwell_coord_temp)));
    end
end

%----------------------------------------------------------------------
if track~='reroute';
    
    referencewell = ones(size(pos,1),1) *-1;
    nonzeros = find(whichTraj>0);
    referencewell(nonzeros) = trajwells(whichTraj(nonzeros),1);
    %assign the times when the animal exited and entered the same well to that
    %reference well (if the time step has not already been defined to a
    %reference well)
    if ~isempty(samewelltraj)
        referencewell(samewelltrajind) = samewelltraj;
    end
    %any remaining -1's in the referencewell vector are times that are
    %too difficult to assign a reference well to
    
    
    %calculate the segment direction for each time step
    segmentdir = ones(size(pos,1),size(lindist,2));
    for i = 1:size(lindist,2)
        segmentdir(:,i) = wellSegmentInfo.segmentDirection(i,segment);
    end
    
    %----------------------------------------------------------------------
    
    % compute head directions relative to the track segments (positive is in outbound direction,
    % negative is inbound)
    vlen = sqrt(vect(:,1).^2 + vect(:,2).^2);
    % normalize vect
    vect(:,1) = vect(:,1) ./ vlen;
    vect(:,2) = vect(:,2) ./ vlen;
    headdir = pos(:,toknum);
    y = sin(headdir);
    x = cos(headdir);
    segheaddir = (vect(:,1) .* x + vect(:,2) .* y);
    segheaddir = repmat(segheaddir,[1 size(lindist,2)]);
    %reverse the head direction for the times when the segment lineaization is
    %reverse relative to the trajectory
    segheaddir(find(segmentdir == 0)) = segheaddir(find(segmentdir == 0))*-1;
    
end
%% Compute velocity
% the following is not good enough since the lindist refers to each
% intersection and change in distance across intersections can't be
% calculated accrately.


% smooth the linear distances with the filter and then go through the linear
% distance positions and take all of the differences to get velocities
% velocity = [];
% for i = 1:size(lindist,2)
%     smoothdist = smoothvect(lindist(:,i), filt);
%     v = diff(smoothdist) / timestep;
%     v = [v(1); v];
%     velocity = [velocity v];
% end

% so do this

% count the difference in linear distance between each time point but reset
% everytime the segment index changes

%index of transitions between segments
Ind(:,1) = [1 find(diff(segment))'+1];
Ind(:,2) = [find(diff(segment))' length(segment)];
 % the default filter for smoothing motion direction is a n second long gaussian
    %with a n/4 second stdev

npoints = smoothwidth/timestep;
filtstd = smoothwidth/(4*timestep);
filt = gaussian(filtstd, npoints);

D=0;
for i=1:size(Ind,1);

    dD= diff(lindist(Ind(i,1):Ind(i,2),1));
    % see if the last transition is close to the start of the reference
    % intersection or away from the reference intersection
    
    if lindist(Ind(i,2))<0.5*segmentInfo.segmentLength(segment(Ind(i,2)));
        lastdistancevalue=lindist(Ind(i,2));
    else
        lastdistancevalue=segmentInfo.segmentLength(segment(Ind(i,2)))-lindist(Ind(i,2));
    end
    
    dD=[dD; lastdistancevalue];
    D=[D; dD];
    
    
end
D=D(2:end,1);

%smoothlindist=smoothvect(D, filt);
velocity=D./timestep;


  

%% collect output

if track~='reroute';
    
    statematrix.time = pos(:,1);
    statematrix.segmentIndex = segment;
    statematrix.wellExitEnter = welltraj;
    statematrix.segmentHeadDirection = segheaddir;
    statematrix.linearVelocity = velocity;
    statematrix.referenceWell = referencewell;
    statematrix.linearDistanceToWells = lindist;
    
    
else
    statematrix.time = pos(:,1);
    statematrix.segmentIndex = segment;
    statematrix.wellExitEnter = welltraj;
    statematrix.segmentHeadDirection = [];
    statematrix.linearVelocity = velocity;
    statematrix.referenceWell = referencewell;
    statematrix.linearDistanceToWells = lindist;
    
end
%%
[statematrix.traj, statematrix.lindist] = gettraj(statematrix, segmenttable, trajwells,track);


function [segmenttable, intersectcoord, wellindex, trajwells, wellSegmentInfo, segmentInfo] = getSegmentTable(firstcoord,track,wellno)
%make a lookup table where each row is [traj trajsegment segID]
%
%well index gives the traj and segment index for each unique trajectory
%endpoint (well location).  If a well is used in multiple trajecories, only
%the index to the first found trajectory is used for that well.
%
%for each trajectory, a row entry in trajwell gives the start and end well
%indeces into wellindex
%
%wellSegmentInfo gives information specific to each well, such as the
%linear distance from that well to each segment.
%
%segmentInfo gives information about each segment, such as segment length
%intersectcoord: add locations of track intersections, as distinct from locations of reward
%wells

tmpcoord = [];
coordtraj = [];
trajlength = [];
wellindex = [];
wells = [];
trajwells = [];




%first, we want to make a list of all the unique wells and all the unique
%segments from the trajectories given
for i = 1:length(firstcoord)
    %add the first and last coordinates of the trajectory as the well
    %locations if these wells have not already been added to the list
    
    if ~rowfind(firstcoord{i}(1,:),wells)
        wellindex = [wellindex; [i 1]];
        wells = [wells; firstcoord{i}(1,:)];
        trajwells(i,1) = size(wellindex,1);
    else
        trajwells(i,1) = rowfind(firstcoord{i}(1,:),wells);
    end
    if ~rowfind(firstcoord{i}(size(firstcoord{i},1),:),wells)
        wellindex = [wellindex; [i size(firstcoord{i},1)]];
        wells = [wells; firstcoord{i}(size(firstcoord{i},1),:)];
        trajwells(i,2) = size(wellindex,1);
    else
        trajwells(i,2) = rowfind(firstcoord{i}(size(firstcoord{i},1),:),wells);
    end
    
    tmpcoord = [tmpcoord ; [firstcoord{i}(1:end-1,:) firstcoord{i}(2:end,:)]]; %get the coordinates of the start and end points of each segment
    tmpcoordtraj = ones(size(firstcoord{i},1)-1,1)*i;
    tmpcoordseg = (1:(size(firstcoord{i},1)-1))';
    trajlength(i,1:2) = [size(firstcoord{i},1)-1 size(coordtraj,1)];
    coordtraj = [coordtraj;[tmpcoordtraj tmpcoordseg]];
end

% intersection coordinates
intersectcoord=wells;

% get wells from Data{day}{epoch}.Wellinfo.rewardedwells
wells=wells(wellno,:);

% get trajwells, the start well and end well, which should be the same for
% all trials of the epoch

trajwells=rowfind(wells(:,1),wells(:,1));
trajwells(:,2)=ones(1,2);

% specify wellindex - override previous loop

wellindex=trajwells;

%any repeating segments are double labeled instead of given a unique
%segment number.
segmentCoords = [];
coordind = [];
indcount = 0;
for i = 1:size(tmpcoord,1)
    [test, testind] = intersect(tmpcoord(i,:),segmentCoords,'rows','legacy');
    reverse = 0;
    if ~(test)
        %if the segment wasn't recorded in one direction, check the reverse
        [test, testind] = intersect([tmpcoord(i,3:4) tmpcoord(i,1:2)],segmentCoords,'rows','legacy');
        reverse = 1;
    end
    origin = tmpcoord(i,1:2);
    endpoint = tmpcoord(i,3:4);
    seglength = sqrt( ((endpoint(1) - origin(1))^2) + ((endpoint(2) - origin(2))^2) );
    if ~(test)
        %this segment hasn't been added to the new list yet, so give it a
        %new index
        segmentCoords = [segmentCoords; tmpcoord(i,:)];
        indcount = indcount+1;
        segmentLength(indcount) = seglength;
        coordind = [coordind; indcount];
        %record which trajectories this segment belongs to
        segmentTrajectories(indcount,1) = bitset(0,coordtraj(i,1),1);
    else
        %this segment is already in the new list, so give it the original
        %index
        coordind = [coordind;testind];
        segmentTrajectories(testind) = bitset(segmentTrajectories(testind),coordtraj(i,1),1);
    end
end
segmentInfo.segmentCoords = segmentCoords;
segmentInfo.segmentLength = segmentLength;
segmentInfo.segmentTraj = segmentTrajectories;
segmenttable = [coordtraj coordind];

%find which segments connect to the start and end of each segment
for i = 1:size(segmentCoords,1)
    tmp = find( ((segmentCoords(:,1) == segmentCoords(i,1))&(segmentCoords(:,2) == segmentCoords(i,2))) | ((segmentCoords(:,3) == segmentCoords(i,1))&(segmentCoords(:,4) == segmentCoords(i,2))) );
    startLinkSegments{i} = setdiff(tmp,i);
    tmp = find( ((segmentCoords(:,1) == segmentCoords(i,3))&(segmentCoords(:,2) == segmentCoords(i,4))) | ((segmentCoords(:,3) == segmentCoords(i,3))&(segmentCoords(:,4) == segmentCoords(i,4))) );
    endLinkSegments{i} = setdiff(tmp,i);
end

%-------------------------------------------------------------------------
if track~='reroute';
    
    %Next we want to create a connectivity matrix stating which segments are
    %directly connected.  This will be used to calculate the distance from any segment to any other segment.
    %Each value in the matrix is defined as e^(-distance) from segA to segB, or
    %e^(-lengthSegA).  This allows us to multiply the values and get summed distance.
    connectivityTable = [];
    wellsegments = [];
    uniqueSegments = unique(segmenttable(:,3))';
    distanceDivisor = 1000;
    if (length(uniqueSegments) > 1)
        for i = uniqueSegments
            segindex = find(segmenttable(:,3) == i)';
            for j = segindex
                traj = segmenttable(j,1);
                segnum = segmenttable(j,2);
                seglength = segmentLength(i);
                %find the segment indeces of the trajectory endpoints
                if (sum((wellindex(:,1) == traj)&((wellindex(:,2) == segnum)|(wellindex(:,2) == (segnum+1)))))
                    transindex = min(find((wellindex(:,1) == traj)&((wellindex(:,2) == segnum)|(wellindex(:,2) == (segnum+1)))));
                    wellsegments = [wellsegments; [i transindex]];
                end
                
                %find the segment indeces of all segments that are directly
                %connected to this segment
                thisSegCoord = segmentCoords(i,:);
                connectedSegments = find( ((segmentCoords(:,1) == thisSegCoord(1))&(segmentCoords(:,2) == thisSegCoord(2))) | ...
                    ((segmentCoords(:,1) == thisSegCoord(3))&(segmentCoords(:,2) == thisSegCoord(4))) | ...
                    ((segmentCoords(:,3) == thisSegCoord(1))&(segmentCoords(:,4) == thisSegCoord(2))) | ...
                    ((segmentCoords(:,3) == thisSegCoord(3))&(segmentCoords(:,4) == thisSegCoord(4))) );
                
                connectedSegments = unique(setdiff(connectedSegments,i));
                %connectedSegments = segmenttable(find((segmenttable(:,1) == traj) & ( (segmenttable(:,2) == segnum+1)|(segmenttable(:,2) == segnum-1) )),3);
                
                distanceDivisor = 1000;
                for k = connectedSegments'
                    %the (i,k)th entry in the connectivity matrix says if the
                    %ith segment is directly connected to the kth segment,
                    %and by what distance: (exp(-lengthOfSegmentI))
                    connectivityTable(i,k) = exp(-seglength/distanceDivisor);
                end
            end
        end
    else
        %there is only one segment
        wellsegments = [1 1];
        connectivityTable = exp(-segmenttable(1,3));
    end
else wellsegments=[1 1];
end

%-------------------------------------------------------------------------

%wellsegments = sortrows(wellsegments,2);
%wellsegments = wellsegments(:,1);
% specify wellsegments based on rewardwell coordinates
wellsegments = rowfind(wells(:,1),segmentCoords(:,3));

%-------------------------------------------------------------------------
%if track~='reroute';

% for i = 1:length(wellsegments)
%     %calculate the length of the arms leading up to the wells (until an
%     %intersection is hit). This will be used later to determine if the
%     %animal completed a trajectory
%     armlength(i) = 0;
%     segcount = 1;
%     foundIntersection = 0;
%     segindex = wellsegments(i);
%     tempTable = connectivityTable;
%     while( (foundIntersection == 0) & (segcount < length(uniqueSegments)) )
%         if(sum(tempTable(segindex,:) > 0) > 1)
%             foundIntersection = 1;
%         else
%             %we have not reached an intersection yet, so find the next
%             %segments in the connection tree
%             tempTable = tempTable*connectivityTable;
%             segcount = segcount + 1;
%         end
%     end
%     %calculate the length of the independant arm (all nonzero numbers on
%     %the row should be the same, so just pick the maximum one)
%     armlength(i) = -log(max(tempTable(segindex,:)));
% end
%end
%------------------------------------------------------------------------

% caculate armlength/distance from well to intersection
armlength=-log(segmentInfo.segmentLength(wellsegments));

%calculate the linear length from each segment to every other segment
%and the sequence of segments to get from segA to segB

if track~='reroute';
    
    tempTable = connectivityTable;
    distanceTable = connectivityTable;
    pathTable = cell(size(connectivityTable));
    diagIndeces = find(eye(size(connectivityTable,1))); %these are the indeces to the diagonal entries
    [nonZerosi,nonZerosj] = find(distanceTable > 0);
    for i = 1:length(nonZerosi)
        pathTable{nonZerosi(i),nonZerosj(i)} = [nonZerosi(i) nonZerosj(i)];
    end
    
    while(sum(distanceTable(:) == 0) > 0)   %we do this loop until we have found a path between every pair of segments
        
        
        oldtempTable = tempTable;
        tempTable = tempTable*connectivityTable;  %this finds the grandchildren...greatgrandchildren...and so on, of each segment, and the corresponding distance
        zerovals = find(distanceTable == 0);
        %which segment pairs became connected?
        [switchedValsi, switchedValsj] = find( (distanceTable == 0) & (tempTable > 0) );
        
        %we only update the values that were zero (because we already found the
        %minimum path for the other pairs)
        distanceTable(zerovals) = tempTable(zerovals);
        
        for i = 1:length(switchedValsi)
            if (switchedValsi(i) ~= switchedValsj(i))
                %find which segment links one path to the new segment
                
                linker = find( (oldtempTable(switchedValsi(i),:) > 0) & (connectivityTable(:,switchedValsj(i))' > 0) );
                prepath = pathTable{switchedValsi(i),linker}(1:end-1);
                postpath = pathTable{linker,switchedValsj(i)}(2:end);
                %add the path to pathTable
                pathTable{switchedValsi(i),switchedValsj(i)} = [prepath linker postpath];
            end
        end
    end
    for i = 1:length(diagIndeces)
        %the diagonal entries are wrong, so we fix them
        distanceTable(diagIndeces(i)) = 1;
        pathTable{diagIndeces(i)} = i;
    end
    
    
    
    %convert the exponent distance back to normal distance
    %distanceTable =(-log(distanceTable)) - log(multtable);
    distanceTable =-log(distanceTable) * distanceDivisor;
    
    %only keep the distances from the wells (we don't care about the other
    %segment pairs)
    distanceTable = distanceTable(wellsegments,:);
    pathTable = pathTable(wellsegments,:);
    segmentDirection = [];
    
    %finally, we need to calculate whether a segment is aligned in the foreward
    %or backward direction relative to each well
    for i = 1:size(pathTable,1)
        for j = 1:size(pathTable,2)
            if ( (isempty(startLinkSegments{j})) & (isempty(endLinkSegments{j})) )
                %there is only one segment
                if (i==1)
                    segmentDirection(i,j) = 1;
                else
                    segmentDirection(i,j) = 0;
                end
            else
                foundpath = 0;
                for k = pathTable{i,j}
                    if (intersect(k,startLinkSegments{j},'legacy'))
                        segmentDirection(i,j) = 1; %from this well, the segment will linearize in the foreward direction
                        foundpath = 1;
                        break;
                    elseif (intersect(k,endLinkSegments{j},'legacy'))
                        segmentDirection(i,j) = 0; %or from the backward direction
                        foundpath = 1;
                        break;
                    end
                end
                if (foundpath == 0)
                    %this is the same segment as the well, and it will
                    %therefore have either no start segments or no end
                    %segments
                    if (isempty(startLinkSegments{j}))
                        segmentDirection(i,j) = 1;
                    elseif (isempty(endLinkSegments{j}))
                        segmentDirection(i,j) = 0;
                    end
                end
            end
        end
    end
    
end

trajwells = [];
for ii = 1:size(segmentInfo.segmentCoords,1)
    trajwells = [trajwells; rowfind(segmentInfo.segmentCoords(ii,1:2), intersectcoord) ...
        rowfind(segmentInfo.segmentCoords(ii,3:4), intersectcoord)];
end

%------------------------------------------------------------------------

%create the wellSegmentInfo structure
if track~='reroute';
    
    wellSegmentInfo.distanceTable = distanceTable;
    wellSegmentInfo.segmentIndex = wellsegments;
    wellSegmentInfo.distanceToIntersection = armlength;
    wellSegmentInfo.wellCoord = wells;
    wellSegmentInfo.pathTable = pathTable;
    wellSegmentInfo.segmentDirection = segmentDirection;
    
else
    wellSegmentInfo.distanceTable = [];
    wellSegmentInfo.segmentIndex = wellsegments;
    wellSegmentInfo.distanceToIntersection = [];
    wellSegmentInfo.wellCoord = wells;
    wellSegmentInfo.pathTable = [];
    wellSegmentInfo.segmentDirection = [];
end
%------------------------------------------------------------
%------------------------------------------------------------


function [trajmatrix trajmatrix_raw] = gettraject(pos, welllocations,trajwells, lindist, wellSegmentInfo, welldist, mindiff, welldistthresh,branchpointdist,track,intersectcoord);
%creates a 3 column matrix, where each row descibes a complete well-to-well
%trajectory. The columns are [starttime startwell endwell]
%welldist is the detection radius around each well (in cm)
%mindiff is the minimum time (in seconds) between detection to be called a trajectory
%mindist is the threshold distance (in cm) for the program to add a trajectory when the animal turns around
%before reaching the end well
%modify: need to get times when rat enters and leaves each junction inorder
%to reconstruct the route used by the rat on each trial

%redefine foodwells as intersection points




nfoodwells = size(intersectcoord,1);
exittimes = [];
minsamples = round(mindiff/(pos(2,1)-pos(1,1)));

% find the times that the animal was within welldist of the foodwells
for i = 1:nfoodwells
    
    tmpwellloc = intersectcoord(i,:);
    % find the times where the animal exited the bounding circle around the well.
    etimes = pos(find(dist(tmpwellloc, pos(:,2:3)) < welldist),1);
    
    if (~isempty(etimes))
        
        % remove all points less than mindiff apart in time
        % this will leave us with one time point per exitpoint indicating when
        % the animal left that exitpoint on each trajectory
        tmpdiff = etimes(2:length(etimes)) - etimes(1:length(etimes) - 1);
        % the leaving times for this exitpoint are the first time of each pair whose
        % difference is larger than mindiff;
        tmp = find(tmpdiff > mindiff);
        
        % tmp omits the last exitpoint, add in the last etime
        tmparray = zeros(length(tmp)+1,2);
        
        tmparray(:,1) = [etimes(tmp) ; etimes(end)];
        tmparray(:,2) = i;
        exittimes = [exittimes ; tmparray];
        
    end
end

if ~isempty(exittimes)
    exittimes = sortrows(exittimes, 1);
end

if track~='reroute';
    
    % this loop detects if the animal went down part more than mindist of an unnoticed arm.
    % if so, then a trajectory is added. Do do this, I use a loop that breaks
    % and starts over every time a new trajecory is added.  This ensures that any other arm
    %entries during that period are also detected.
    
    foundone = 0;
    stopcue = 0;
    
    
    while(stopcue == 0)
        foundone2 = 0;
        for ct = 2:size(exittimes,1)
            tmptraj = [exittimes(ct-1,2) exittimes(ct,2)];
            tmptimes = [exittimes(ct-1,1) exittimes(ct,1)];
            tmpindex = find((pos(:,1)> tmptimes(1)) & (pos(:,1) < tmptimes(2)));
            if (~isempty(tmpindex))
                for checkout = 1:size(lindist,2)
                    if ~intersect(checkout, tmptraj,'legacy')
                        if isfinite(branchpointdist)
                            %use the distance from the branch point
                            inarmind = find( (lindist(tmpindex,checkout) < (wellSegmentInfo.distanceToIntersection(checkout))-branchpointdist));
                        else
                            %use the threshhold distance from the well
                            inarmind = find( (lindist(tmpindex,checkout) < (wellSegmentInfo.distanceToIntersection(checkout))) & ...
                                (lindist(tmpindex,checkout) < welldistthresh) );
                        end
                        if ~isempty(inarmind)
                            if (length(inarmind) > 1)
                                [tmpmin,endind] = min(find((diff(inarmind)>=minsamples)));
                                if isempty(endind)
                                    endind = length(inarmind);
                                end
                            else
                                endind = 1;
                            end
                            [minval, minind] = min(lindist(tmpindex(inarmind(1:endind)),checkout));
                            outtime = pos(tmpindex(inarmind(minind)),1);
                            exittimes = [exittimes(1:ct-1,:);outtime checkout;exittimes(ct:end,:)];
                            foundone = 1;
                            foundone2 = 1;
                            break;
                        end
                    end
                end
            end
            if (foundone)
                %a trajectory was added, so we break and start over
                foundone = 0;
                break;
            end
        end
        if (foundone2 == 0)
            %no trajectory was added, so we are done
            stopcue = 1;
        else
            foundone2 = 0;
        end
    end
    
end
%create a matrix describing the start time, start well, and end well of
%each complete trajectory
trajmatrix = [];
for i = 1:size(exittimes,1)-1
    trajmatrix(i,1:3) = [exittimes(i,1) exittimes(i,2) exittimes(i+1,2)];
end

% find segmentation errors, where an intersection is skipped
% this means the rat's trajectory is greater than welldist from one of the
% intersections the intersection was not calculated as crossed

testmatrix=sum(trajmatrix(2:end,2)-trajmatrix(1:end-1,3));
if testmatrix~=0
    disp('need to redo task structure and ensure intersection zones contain trajectories');
end




trajmatrix_raw = trajmatrix;

%removing duplicate well entries
ind_to_remove = find((trajmatrix(:,2)-trajmatrix(:,3))==0);
trajmatrix(ind_to_remove,:)=[];










function [state, lindist] = gettraj(statematrix, segmenttable, trajwells,track)

%get the most probable linear trajectory number, based on well exit/enters
%and head direction.

if track~='reroute';
    
    traject = trajwells;
    includeStates = 1:6;
    
    %if no reference well is assigned, it is not a valid point
    %also, only points inside the designated time range are valid
    validPoints = ((statematrix.referenceWell > 0));
    
    validPointsIndex = find(validPoints);
    trajvector = ones(length(validPointsIndex),1)*-1;
    statevector = ones(length(validPointsIndex),1)*-1;
    referencewell = statematrix.referenceWell(validPointsIndex);
    welltraj = statematrix.wellExitEnter(validPointsIndex,:); %exit and enter wells
    segmentindex = statematrix.segmentIndex(validPointsIndex);
    samewell = (welltraj(:,1) == welltraj(:,2)); %Which times are during trajectories between the same wells
    diffwell = (welltraj(:,1) ~= welltraj(:,2)); %Which times are during trajectories between different wells
    trajcount = rowcount(statematrix.segmentIndex(validPointsIndex),segmenttable(:,3)); %how many of the linear trajectories does each segment belong to
    
    matrixIndex = sub2ind(size(statematrix.linearDistanceToWells),validPointsIndex,referencewell);
    lineardistance = statematrix.linearDistanceToWells(matrixIndex);
    
    forewarddir = ((statematrix.segmentHeadDirection(matrixIndex) >= 0)); %facing positve dir and moving positive dir
    backwarddir = ((statematrix.segmentHeadDirection(matrixIndex) < 0));%facing negative dir and moving negative dir
    
    
    
    
    %LEVEL 1
    for trajnum = 1:size(traject,1)
        currtraj = traject(trajnum,:);
        
        %which times occur when the animal is intransit between the two wells
        %for the current trajectory
        inDefinedTraj(:,trajnum) = (((welltraj(:,1)==currtraj(1))&(welltraj(:,2)==currtraj(2))) | ((welltraj(:,1)==currtraj(2))&(welltraj(:,2)==currtraj(1))));
        %which segments are in the current trajectory
        segmentsInTraj = segmenttable(find(segmenttable(:,1) == trajnum),3);
        
        %find all times when the animal was on the current trajectory (either
        %in the foreward or backward direction).  This is defined as the times
        %when the animal is in transit between the two correct endpoints and is
        %on a track segment that is part of the trajectory.
        
        %first find the foreward direction times (both head and motion
        %direction in the positive direction, and velocity greater than
        %minvelocity)
        findindex = find( inDefinedTraj(:,trajnum)& intersect(segmentindex,segmentsInTraj,'legacy') & forewarddir);
        %assign the current trajectory to those indeces
        trajvector(findindex) = 2*(trajnum)-1;
        if intersect(1,includeStates,'legacy')
            statevector(findindex) = trajvector(findindex);
        end
        %then find the negative direction times
        findindex = find( inDefinedTraj(:,trajnum) & intersect(segmentindex,segmentsInTraj,'legacy') & backwarddir);
        trajvector(findindex) = 2*(trajnum);
        if intersect(1,includeStates,'legacy')
            statevector(findindex) = trajvector(findindex);
        end
    end
    
    undefinedindex = (trajvector == -1); %still undefined times
    inAnyDefinedTraj = (sum(inDefinedTraj,2) > 0); %These times occur when the animal is in transit between two wells of a defined trajectory
    
    % LEVEL 2
    % for the times that are still undefined and not occurring during any of the above well-to-well trajectories,
    % and occur when the animal is not leaving and entering the same well,
    % and is in a track segment that only belongs to one of the defined linear
    % trajectories, assign the proper trajectory
    findindex = find((undefinedindex) & (diffwell) & (trajcount == 1) & (forewarddir) & (~inAnyDefinedTraj));
    trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1))-1;
    if intersect(2,includeStates,'legacy')
        statevector(findindex) = trajvector(findindex);
    end
    findindex = find((undefinedindex) & (diffwell) & (trajcount == 1) & (backwarddir) & (~inAnyDefinedTraj));
    trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1));
    if intersect(2,includeStates,'legacy')
        statevector(findindex) = trajvector(findindex);
    end
    
    if ((nargout == 1) & (max(includeStates) <= 2))
        return
    end
    
    % LEVEL 3
    % fill all times when the animal is on a segment that belongs to only one
    % trajectory (even if it is leaving and entering the same well)
    undefinedindex = (trajvector == -1); %still undefined times
    findindex = find((undefinedindex) & (trajcount == 1) & (forewarddir));
    trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1))-1;
    if intersect(3,includeStates,'legacy')
        statevector(findindex) = trajvector(findindex);
    end
    findindex = find((undefinedindex) & (trajcount == 1) & (backwarddir));
    trajvector(findindex) = 2*(segmenttable(rowfind(segmentindex(findindex),segmenttable(:,3)),1));
    if intersect(3,includeStates,'legacy')
        statevector(findindex) = trajvector(findindex);
    end
    
    % LEVEL 4
    % Fill in the undefined times when the animal was on an amiguous segment (like the
    % home arm in a w-maze)
    undefinedindex = (trajvector == -1); %the undefined times
    %if facing positive direction, then assign the traj of the next traj in the
    %future
    findindex = find((undefinedindex) & (trajcount > 1) & (forewarddir));
    trajvector = vectorfill(trajvector, -1, 1, findindex);
    findex2 = findindex(find(~mod(trajvector(findindex),2)));
    trajvector(findex2) = trajvector(findex2)-1;
    if intersect(4,includeStates,'legacy')
        statevector(findindex) = trajvector(findindex);
    end
    %otherwise assign the closest traj that happened in the past
    findindex = find((undefinedindex) & (trajcount > 1) & (backwarddir));
    trajvector = vectorfill(trajvector, -1, -1, findindex);
    findex2 = findindex(find(mod(trajvector(findindex),2)));
    trajvector(findex2) = trajvector(findex2)+1;
    if intersect(4,includeStates,'legacy')
        statevector(findindex) = trajvector(findindex);
    end
    
    % LEVEL 5
    % now we may have some undefined points in the beginning and end of the
    % session- we assign these to the first trajectory
    undefinedindex = (trajvector < 0); %still undefined times
    trajvector(find(undefinedindex)) = -1;
    statevector(find(undefinedindex)) = -1;
    findindex = find((undefinedindex) & (forewarddir));
    trajvector(findindex) = 1;
    if intersect(5,includeStates,'legacy')
        statevector(findindex) = trajvector(findindex);
    end
    findindex = find((undefinedindex) & (backwarddir));
    trajvector(findindex) = 2;
    if intersect(5,includeStates,'legacy')
        statevector(findindex) = trajvector(findindex);
    end
    
    % LEVEL 6
    % Fill in all remaining points, regardless of running velocity
    undefinedindex = find(trajvector == -1); %still undefined times
    trajvector = vectorfill(trajvector, -1, 0, undefinedindex);
    if intersect(6,includeStates,'legacy')
        statevector = trajvector;
    end
    
    
    
    
    state = ones(length(statematrix.time),1)*-1;
    lindist = ones(length(statematrix.time),1)*-1;
    state(validPointsIndex) = statevector;
    lindist(validPointsIndex) = lineardistance;
    
else
    
    state = ones(length(statematrix.time),1)*-1;
    lindist = [];
    %state(validPointsIndex) = statevector;
    %lindist(validPointsIndex) = lineardistance;
end




