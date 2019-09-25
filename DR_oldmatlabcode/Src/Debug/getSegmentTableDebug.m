
function [segmenttable, wellindex, trajwells, wellSegmentInfo, segmentInfo] = getSegmentTableDebug(firstcoord,track)
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

tmpcoord = [];
coordtraj = [];
trajlength = [];
wellindex = [];
wells = [];
trajwells = [];

%first, we want to make a list of all the unique wells and all the unique
%segments from the trajectories given
for i = 1:length(firstcoord) % firstcoord contains the start end end points of each segment [x y]
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
         trajwells(i,2) = rowfind(firstcoord{i}(size(firstcoord{i},1),:),wells); % start and end points of a trajectory based on intersection points
     end
          
     tmpcoord = [tmpcoord ; [firstcoord{i}(1:end-1,:) firstcoord{i}(2:end,:)]]; %get the coordinates of the start and end points of each segment
     tmpcoordtraj = ones(size(firstcoord{i},1)-1,1)*i;
     tmpcoordseg = (1:(size(firstcoord{i},1)-1))';
     trajlength(i,1:2) = [size(firstcoord{i},1)-1 size(coordtraj,1)];
     coordtraj = [coordtraj;[tmpcoordtraj tmpcoordseg]];
end

% respecify wells as four end arms of the track (points 6 7 8 9) for reroute task
if strcmp(track,'reroute')==1;
    wells=wells([6 7 8 9],:);
end

%any repeating segments are double labeled instead of given a unique
%segment number.
segmentCoords = [];
coordind = [];
indcount = 0;
for i = 1:size(tmpcoord,1)
    [test, testind] = ismember(tmpcoord(i,:),segmentCoords,'rows');
    reverse = 0;
    if ~(test)
        %if the segment wasn't recorded in one direction, check the reverse
        [test, testind] = ismember([tmpcoord(i,3:4) tmpcoord(i,1:2)],segmentCoords,'rows');
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
segmentInfo.segmentCoords = segmentCoords; % coordinate of each segment [x1 y1 x2 y2]
segmentInfo.segmentLength = segmentLength; % length of each segment
segmentInfo.segmentTraj = segmentTrajectories; %?
segmenttable = [coordtraj coordind];

%find which segments connect to the start and end of each segment
for i = 1:size(segmentCoords,1)
    tmp = find( ((segmentCoords(:,1) == segmentCoords(i,1))&(segmentCoords(:,2) == segmentCoords(i,2))) | ((segmentCoords(:,3) == segmentCoords(i,1))&(segmentCoords(:,4) == segmentCoords(i,2))) );
    startLinkSegments{i} = setdiff(tmp,i);   
    tmp = find( ((segmentCoords(:,1) == segmentCoords(i,3))&(segmentCoords(:,2) == segmentCoords(i,4))) | ((segmentCoords(:,3) == segmentCoords(i,3))&(segmentCoords(:,4) == segmentCoords(i,4))) );
    endLinkSegments{i} = setdiff(tmp,i); % each column number refers to one segment, contents of that column refers to all connected segments
end

%-------------------------------------------------------------------------
%if track~='reroute';

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
%else wellsegments=[1 1];
%end

%-------------------------------------------------------------------------

wellsegments = sortrows(wellsegments,2);
wellsegments = wellsegments(:,1);

%-------------------------------------------------------------------------
%if track~='reroute';

for i = 1:length(wellsegments)
    %calculate the length of the arms leading up to the wells (until an
    %intersection is hit). This will be used later to determine if the
    %animal completed a trajectory
    armlength(i) = 0;
    segcount = 1;
    foundIntersection = 0;
    segindex = wellsegments(i);
    tempTable = connectivityTable;
    while( (foundIntersection == 0) & (segcount < length(uniqueSegments)) )
        if(sum(tempTable(segindex,:) > 0) > 1)
            foundIntersection = 1;            
        else
            %we have not reached an intersection yet, so find the next
            %segments in the connection tree
            tempTable = tempTable*connectivityTable;
            segcount = segcount + 1;
        end
    end
    %calculate the length of the independant arm (all nonzero numbers on
    %the row should be the same, so just pick the maximum one)           
    armlength(i) = -log(max(tempTable(segindex,:)));      
end
%end
%------------------------------------------------------------------------


%calculate the linear length from each segment to every other segment
%and the sequence of segments to get from segA to segB

%if track~='reroute';

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
                if (ismember(k,startLinkSegments{j}))
                    segmentDirection(i,j) = 1; %from this well, the segment will linearize in the foreward direction
                    foundpath = 1;
                    break;
                elseif (ismember(k,endLinkSegments{j}))
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

%end

%------------------------------------------------------------------------

%create the wellSegmentInfo structure
%if track~='reroute';
    
wellSegmentInfo.distanceTable = distanceTable;
wellSegmentInfo.segmentIndex = wellsegments;
wellSegmentInfo.distanceToIntersection = armlength; 
wellSegmentInfo.wellCoord = wells;
wellSegmentInfo.pathTable = pathTable;
wellSegmentInfo.segmentDirection = segmentDirection;   

%else 
wellSegmentInfo.distanceTable = [];
wellSegmentInfo.segmentIndex = wellsegments;
wellSegmentInfo.distanceToIntersection = []; 
wellSegmentInfo.wellCoord = wells;
wellSegmentInfo.pathTable = [];
wellSegmentInfo.segmentDirection = [];   
%end

end

