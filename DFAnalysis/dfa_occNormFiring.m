function [out] = dfa_occNormFiring(index, excludeperiods, spikes,linpos, pos, ...
    task, varargin)
%
%Calculates the 2d occupancy normalized firing rate for the cell and
%organizes the output into the different trajectories.
%
%binsize- the length of each spatial bin (default 1 cm)
%std - defines the shape of the 2d gaussian used to smooth spikerate.
%              (default 1)
%The output is a cell array where each cell contains a cell
%descibing one trajectory


% parse the options
appendindex = 1;
binsize = 1; % cm squarec
stdocc = 3;
stdspike = 3;
threshocc = 0.002; % Threshold occupancy in seconds

if (~isempty(varargin))
    assign(varargin{:});
end

warning('OFF','MATLAB:divideByZero');
out.spikes = [];
out.spikerate = [];
out.smoothedspikerate = [];
out.occupancy = [];
out.smoothedoccupancy = [];
out.xticks = [];
out.yticks = [];
out.seqX = [];
out.allseqX = [];
out.allseqY = [];
out.seqY = [];
out.allseqX1 = [];
out.allseqY1 = [];
out.wellCoords = [];
out.index = index;
try spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
catch
    return
end

spikesData = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
posdata = pos{index(1)}{index(2)}.data;
posfields = pos{index(1)}{index(2)}.fields;
taskEnv = task{index(1)}{index(2)}.environment;

if isempty(spikesData)
    fprintf('spikes empty \n');
    return
end    
    

timestring = 'time';
timecol = find(cell2mat(cellfun(@(x) strcmp(x,timestring), strsplit(posfields,...
    ' '), 'UniformOutput', false)));
postime = posdata(:,timecol);

if strcmp(taskEnv, 'wtrack')
    statematrix = linpos{index(1)}{index(2)}.statematrix;
else %if strcmp(taskEnv, 'openfield')
    statematrix.time = postime;
    statematrix.traj = ones(length(postime),1);
end

statevector = statematrix.traj;

% front padding the posdata vect to match statevec
if length(statevector(:,1)) ~= length(posdata(:,1));
    posdata = [zeros((length(statevector(:,1))-length(posdata(:,1))), ...
        length(posdata(1,:))); posdata];
end

% Use Exclude Periods for TimeFilter version in addition to statevector=-1
% Based on exclude time,
statevector(find(isExcluded(posdata(:,1), excludeperiods))) = -1;


posindexfield = knnsearch(postime, spikesData(:,1));

xstring = 'x-loess';
xcol = find(cell2mat(cellfun(@(x) strcmp(x,xstring), strsplit(posfields, ' '), ...
    'UniformOutput', false)));
posxdata = posdata(posindexfield,xcol);
ystring = 'y-loess';
ycol = find(cell2mat(cellfun(@(x) strcmp(x,ystring), strsplit(posfields, ' '), ...
    'UniformOutput', false)));
posydata = posdata(posindexfield,ycol);

if ~isempty(spikesData)
    spikesData = [spikesData(:,1) posxdata posydata posindexfield];
    spikesData(:,5) = statevector(posindexfield); %add the traj state for each spike
else
    spikesData = [0 0 -1];
end

trajnum = max(statevector);
timestep = statematrix.time(2,1) - statematrix.time(1,1);
goodspikes = [];
goodspikeind = (spikesData(:,5) ~= -1);
%create a list of the non sharp-wave spikes
goodspikes = spikesData(goodspikeind,:);
%make a cell array, where each cell contains data for one trajectory.
trajdata = cell(1,trajnum);
goodlocationind = (find(statevector ~= -1 ));
goodlocations = [statematrix.time(goodlocationind) posdata(goodlocationind,[2 3]) statevector(goodlocationind)]; %CHANGED the positionindex at all valid times

out.allseqX(:,1) = 1:(((ceil(max(goodlocations(:,2)))-floor(min(goodlocations(:,2))))/binsize)+1);
out.allseqY(:,1) = 1:(((ceil(max(goodlocations(:,3)))-floor(min(goodlocations(:,3))))/binsize)+1);
out.allseqX(:,2) = floor(min(goodlocations(:,2))):binsize:ceil(max(goodlocations(:,2)));
out.allseqY(:,2) = floor(min(goodlocations(:,3))):binsize:ceil(max(goodlocations(:,3)));
out.allseqX1(:,2) = 1:binsize:ceil(max(goodlocations(:,2)));
out.allseqY1(:,2) = 1:binsize:ceil(max(goodlocations(:,3)));

for i = 1:length(trajdata)
    
    %get all the position indices and spikes when the animal was on the ith
    %trajectory
    tmpposition = goodlocations(find(goodlocations(:,4) == i),[2 3]);
    tmpspikes = goodspikes(find(goodspikes(:,5) == i),[2 3]);
    
    if ~isempty(tmpposition)
        minx = floor(min(tmpposition(:,1)));
        maxx = ceil(max(tmpposition(:,1)));
        binx = (minx:binsize:maxx);
        miny = floor(min(tmpposition(:,2)));
        maxy = ceil(max(tmpposition(:,2)));
        biny = (miny:binsize:maxy);
     
%do this instead of calling hist2.. much faster
        xr = interp1(binx, 1:numel(binx), tmpposition(:,1), 'nearest');
        yr = interp1(biny, 1:numel(biny), tmpposition(:,2), 'nearest');
        xryr = [xr yr];
        xryr = xryr(~any(isnan(xryr),2),:);
        out.occupancy{i} = accumarray( [xryr] , 1, [numel(binx) numel(biny)]);
        out.xticks{i} = binx;
        out.yticks{i} = biny;
        
        
        % create list as long as ticks for anchoring/plotting with a bin size more than 1
        out.seqX{i} = [out.allseqX(find(out.xticks{i}(1) == out.allseqX(:,2)),1):length(out.xticks{i})];
        out.seqY{i} = [out.allseqY(find(out.yticks{i}(1) == out.allseqY(:,2)),1):length(out.yticks{i})];
                
        if ~isempty(tmpspikes)
        
        xr = interp1(binx, 1:numel(binx), tmpspikes(:,1), 'nearest');
        yr = interp1(biny, 1:numel(biny), tmpspikes(:,2), 'nearest');
        xryr = [xr yr];
        xryr = xryr(~any(isnan(xryr),2),:);
        out.spikes{i} = accumarray( [xryr] , 1, [numel(binx) numel(biny)]);
            
            nonzero = find(out.occupancy{i} ~= 0);
            out.spikerate{i} = zeros(size(out.spikes{i}));
            out.spikerate{i}(nonzero) = out.spikes{i}(nonzero) ./(timestep* out.occupancy{i}(nonzero) );
            
             %2) Smooth spikerate and occupancy
             gspike = gaussian2(stdspike,2*stdspike);
            out.smoothedspikerate{i} = filter2(gspike,(out.spikerate{i}));

            gocc = gaussian2(stdocc,2*stdocc);
            out.smoothedoccupancy{i} = filter2(gocc, out.occupancy{i});
            
            % 3) Turn occupancy to seconds and set spikerate wherever occupancy
            %is < threshold occupancy in seconds to 0
            
            out.occupancy{i} = timestep*out.occupancy{i};
            out.smoothedoccupancy{i} = timestep*out.smoothedoccupancy{i};
            
            if any(out.smoothedoccupancy{i} < 0)
                pause
            end
            
            %zero = find(smoothedoccupancy == 0);
            zero = find(out.smoothedoccupancy{i} <= threshocc);
            %output.smoothedspikerate{i}(zero) = -2; %no occupancy is negative and diff/darker from occupancy but no spikes
            out.smoothedspikerate{i}(zero) = -1;   
            
            %return coordinates (from coordprogram) for the well enter exit for this traj
            % right now this is a hack to get w track working
            if strcmp(taskEnv, 'wtrack')
                wellCoords = linpos{index(1)}{index(2)}.wellSegmentInfo.wellCoord;
                traj2wells = [1 1 2; 2 2 1; 3 1 3; 4 3 1]; % traj enterwellID exitwellID
                out.wellCoords{i} = wellCoords(traj2wells(i,2:3),:); % to do [enterX enterY ; exitX exitY]
            else
                out.wellCoords{i} = [];
            end
        else
            out.spikes{i} = [];
            out.spikerate{i} = [];
            out.smoothedspikerate{i} = [];
            out.smoothedoccupancy{i} = [];
            out.wellCoords{i} = [];
        end
        
    else
        out.wellCoords{i} = [];
        out.spikes{i} = [];
        out.spikerate{i} = [];
        out.smoothedspikerate{i} = [];
        out.occupancy{i} = [];
        out.smoothedoccupancy{i} = [];
        out.xticks{i}  = [];
        out.yticks{i} = [];
        out.seqX{i} = [];
        out.allseqX{i} = [];
        out.allseqY{i} = [];
        out.seqY{i} = [];
        out.allseqX1{i} = [];
        out.allseqY1{i} = [];
        
    end

end
if appendindex == 1
        out.index = index;
end

warning('ON','MATLAB:divideByZero');


