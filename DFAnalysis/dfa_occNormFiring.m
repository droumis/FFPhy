function [out] = dfa_occNormFiring(index, excludeperiods, spikes,linpos, pos, task, varargin)
%[output] = twodoccupancy(spikes,statevector, linpos, pos, index, binsize)
%[output] = twodoccupancy(spikes,statevector, linpos, pos, index)
%
%Calculates the 2d occupancy normalized firing rate for the cell and
%organizes the output into the different trajectories.
%
% Separates by traj is the difference with DFAsj_openfieldrate_tf
% Ther is also a twodoccupancy in FUnctions
%
%spikes - the 'spikes' cell array for the day you are analyzing
%linpos - the output of LINEARDAYPROCESS for the day you are analyzing.
%pos - the output of nspike_fixpos
%index - [day epoch tetrode cell]
%statevector - the outputs of GETBEHAVESTATE. This is a vector with the traj
%              number for each position time (1 based). -1 values signify
%              invalid times and are not used.
%binsize- the length of each spatial bin (default 1 cm)
%std - defines the shape of the 2d gaussian used to smooth spikerate.
%              (default 1)
%
%The output is a cell array where each cell contains a cell
%descibing one trajectory.  These cells contain n matrices. The matrices are: occupancy,
%bin vector x, bin vector y, bin spike count, occ normailized firing per bin, and smoothed
%occ normalized firing. If the cell is empty, the animal did not enter that trajectory.
%


% parse the options
appendindex = 1;
binsize = 1; % cm squarec
stdocc = 3;
stdspike = 3;
threshocc = 0.002; % Threshold occupancy in seconds

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
            case 'std'
                std = varargin{option+1};
            case 'threshocc'
                std = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end


warning('OFF','MATLAB:divideByZero');


% segmentInfo = linpos{index(1)}{index(2)}.segmentInfo;
try spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
catch
    output.spikes = [];
    output.spikerate = [];
    output.smoothedspikerate = [];
    output.occupancy = [];
    output.smoothedoccupancy = [];
    output.xticks = [];
    output.yticks = [];
    output.seqX = [];
    output.allseqX = [];
    output.allseqY = [];
    output.seqY = [];
    output.allseqX1 = [];
    output.allseqY1 = [];
    output.wellCoords = [];
    output.index = index;
    out = output;
    return
end



spikesData = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
posdata = pos{index(1)}{index(2)}.data;
posfields = pos{index(1)}{index(2)}.fields;
taskEnv = task{index(1)}{index(2)}.environment;

timestring = 'time';
timecol = find(cell2mat(cellfun(@(x) strcmp(x,timestring), strsplit(posfields, ' '), 'UniformOutput', false)));
postime = posdata(:,timecol);

if strcmp(taskEnv, 'wtrack')
    statematrix = linpos{index(1)}{index(2)}.statematrix;
elseif strcmp(taskEnv, 'openfield')
    statematrix.time = postime;
    statematrix.traj = ones(length(postime),1);
end

statevector = statematrix.traj;

%DR added 4/22/14... front padding the posdata vect to match statevec
if length(statevector(:,1)) ~= length(posdata(:,1));
    posdata = [zeros((length(statevector(:,1))-length(posdata(:,1))), length(posdata(1,:))); posdata];
end

% Use Exclude Periods for TimeFilter version in addition to statevector=-1
statevector(find(isExcluded(posdata(:,1), excludeperiods))) = -1; % Based on exclude time,




posindexfield = knnsearch(postime, spikesData(:,1));

xstring = 'x-sm';
xcol = find(cell2mat(cellfun(@(x) strcmp(x,xstring), strsplit(posfields, ' '), 'UniformOutput', false)));
posxdata = posdata(posindexfield,xcol);
ystring = 'y-sm';
ycol = find(cell2mat(cellfun(@(x) strcmp(x,ystring), strsplit(posfields, ' '), 'UniformOutput', false)));
posydata = posdata(posindexfield,ycol);

% spikes.data fields: time x y dir not_used amplitude(highest variance channel) posindex x-sm y-sm dir-sm
% if size(spikes,2)>7
%     posx = 8; posy = 9;
% else
%     posx = 2; posy = 3;
% end
% posx = 2; posy = 3;
% posindexfield = 7;

if ~isempty(spikesData)
    spikesData = [spikesData(:,1) posxdata posydata posindexfield];
    spikesData(:,5) = statevector(posindexfield); %add the traj state for each spike
    %posindexfield = isdatafield(spikesfields,'posindex');
    
else
    spikesData = [0 0 -1];
end

trajnum = max(statevector);
timestep = statematrix.time(2,1) - statematrix.time(1,1);
goodspikes = [];
goodspikeind = (spikesData(:,5) ~= -1);
%create a list of the non sharp-wave spikes [time traj linearloc]
goodspikes = spikesData(goodspikeind,:);
%make a cell array, where each cell contains data for one trajectory.
%inside each cell [binlocation occupancy spikecount firingrate]

trajdata = cell(1,trajnum);
goodlocationind = (find(statevector ~= -1 ));
goodlocations = [statematrix.time(goodlocationind) posdata(goodlocationind,[2 3]) statevector(goodlocationind)]; %CHANGED the positionindex at all valid times

%     if isempty(goodspikes) %there were no useable spikes
%         goodspikes = [];
%         trajdata = [];
%     end

% Shantanu: out has to be a struct - so I have to change how data is stored
% for all trajectories

%DR
output.allseqX(:,1) = 1:(((ceil(max(goodlocations(:,2)))-floor(min(goodlocations(:,2))))/binsize)+1);
output.allseqY(:,1) = 1:(((ceil(max(goodlocations(:,3)))-floor(min(goodlocations(:,3))))/binsize)+1);
output.allseqX(:,2) = floor(min(goodlocations(:,2))):binsize:ceil(max(goodlocations(:,2)));
output.allseqY(:,2) = floor(min(goodlocations(:,3))):binsize:ceil(max(goodlocations(:,3)));
output.allseqX1(:,2) = 1:binsize:ceil(max(goodlocations(:,2)));
output.allseqY1(:,2) = 1:binsize:ceil(max(goodlocations(:,3)));

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
         % Shantanu - I am using hist2 in ~/Src/Matlab/sj_Utility
%         [output.occupancy{i} output.xticks{i} output.yticks{i}] = hist2(tmpposition(:,1), tmpposition(:,2), binx, biny);
        

%DR do this instead of calling hist2.. this is much faster
        xr = interp1(binx, 1:numel(binx), tmpposition(:,1), 'nearest');
        yr = interp1(biny, 1:numel(biny), tmpposition(:,2), 'nearest');
        xryr = [xr yr];
        xryr = xryr(~any(isnan(xryr),2),:);
        output.occupancy{i} = accumarray( [xryr] , 1, [numel(binx) numel(biny)]);
        output.xticks{i} = binx;
        output.yticks{i} = biny;
        
        
                %DR sometimes hist2 messed up and addes an extra non-integer edge at the end of x edges... recrop matrix according to binx, biny
%         output.occupancy{i} = output.occupancy{i}( 1:length(biny), 1:length(binx));
%         output.xticks{i}  = output.xticks{i}(1,1:length(binx));
%         output.yticks{i} = output.yticks{i}(1,1:length(biny));
        
        %DR create list as long as ticks for anchoring/plotting with a bin size more than 1
        output.seqX{i} = [output.allseqX(find(output.xticks{i}(1) == output.allseqX(:,2)),1):length(output.xticks{i})];
        output.seqY{i} = [output.allseqY(find(output.yticks{i}(1) == output.allseqY(:,2)),1):length(output.yticks{i})];
        
        
        if ~isempty(tmpspikes)
            
             %1) Get spikerate from spikes and occupancy
%             [output.spikes{i} BX BY] = hist2(tmpspikes(:,1), tmpspikes(:,2), binx, biny);
            
           %DR do this instead of calling hist2.. this is much faster
        xr = interp1(binx, 1:numel(binx), tmpspikes(:,1), 'nearest');
        yr = interp1(biny, 1:numel(biny), tmpspikes(:,2), 'nearest');
        xryr = [xr yr];
        xryr = xryr(~any(isnan(xryr),2),:);
        output.spikes{i} = accumarray( [xryr] , 1, [numel(binx) numel(biny)]);
            
                        %DR sometimes hist2 messed up and adds an extra non-integer edge at the end of edges... recrop matrix according to binx, biny
%             output.spikes{i} = output.spikes{i}(1:length(biny), 1:length(binx));
%             BX  = BX(1,1:length(binx));
%             BY = BY(1,1:length(biny));
            
            nonzero = find(output.occupancy{i} ~= 0);
            output.spikerate{i} = zeros(size(output.spikes{i}));
            output.spikerate{i}(nonzero) = output.spikes{i}(nonzero) ./(timestep* output.occupancy{i}(nonzero) );
            
             %2) Smooth spikerate and occupancy

            
             gspike = gaussian2(stdspike,2*stdspike);
            output.smoothedspikerate{i} = filter2(gspike,(output.spikerate{i}));
            %smoothedoccupancy = [];
            %smoothedoccupancy = zeros(size(output(i).spikes));
            gocc = gaussian2(stdocc,2*stdocc);
            output.smoothedoccupancy{i} = filter2(gocc, output.occupancy{i});
            
            % 3) Turn occupancy to seconds and set spikerate wherever occupancy
            %is < threshold occupancy in seconds to 0
            
            output.occupancy{i} = timestep*output.occupancy{i};
            output.smoothedoccupancy{i} = timestep*output.smoothedoccupancy{i};
            
            %zero = find(smoothedoccupancy == 0);
            zero = find(output.smoothedoccupancy{i} <= threshocc);
            %output.smoothedspikerate{i}(zero) = -2; %no occupancy is negative and diff/darker from occupancy but no spikes
            output.smoothedspikerate{i}(zero) = -1;   
            
            %return coordinates (from coordprogram) for the well enter exit for this traj
            % right now this is a hack to get w track working
            if strcmp(taskEnv, 'wtrack')
                wellCoords = linpos{index(1)}{index(2)}.wellSegmentInfo.wellCoord;
                traj2wells = [1 1 2; 2 2 1; 3 1 3; 4 3 1]; % traj enterwellID exitwellID
                output.wellCoords{i} = wellCoords(traj2wells(i,2:3),:); % to do [enterX enterY ; exitX exitY]
            else
                output.wellCoords{i} = [];
            end
        else
            output.spikes{i} = [];
            output.spikerate{i} = [];
            output.smoothedspikerate{i} = [];
            output.smoothedoccupancy{i} = [];
            output.wellCoords{i} = [];
        end
        
    else
        output.wellCoords{i} = [];
        output.spikes{i} = [];
        output.spikerate{i} = [];
        output.smoothedspikerate{i} = [];
        output.occupancy{i} = [];
        output.smoothedoccupancy{i} = [];
        output.xticks{i}  = [];
        output.yticks{i} = [];
        output.seqX{i} = [];
        output.allseqX{i} = [];
        output.allseqY{i} = [];
        output.seqY{i} = [];
        output.allseqX1{i} = [];
        output.allseqY1{i} = [];
        
    end

end
if appendindex == 1
        output.index = index;
end

out = output; 
% Make it a 1X1 cell - Access out{ep}.smoothedspikerate{traj}, for trajs 1,2,3,4
% where ep = epochfilter - usually 1

% if appendindex == 0
%     out = output;
% elseif appendindex == 1
%     output.index = index;
%     out = output;
% end

warning('ON','MATLAB:divideByZero');


