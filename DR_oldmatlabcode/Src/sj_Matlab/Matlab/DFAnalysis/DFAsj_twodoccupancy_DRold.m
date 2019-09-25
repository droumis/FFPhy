function [out] = DFAsj_twodoccupancy(index, excludeperiods, spikes,linpos, pos, varargin)
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
appendindex = 2;
binsize = 2; % cm square
std = 1;
threshocc = 0.02; % Threshold occupancy in seconds

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

% Linpos, spikes and pos
statematrix = linpos{index(1)}{index(2)}.statematrix;
spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
spikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
posdata = pos{index(1)}{index(2)}.data;
statevector = statematrix.traj;

%DR added 4/22/14... front padding the posdata vect to match statevec
if length(statevector(:,1)) ~= length(posdata(:,1));
    posdata = [zeros((length(statevector(:,1))-length(posdata(:,1))), length(posdata(1,:))); posdata];
end


% Use Exclude Periods for TimeFilter version in addition to statevector=-1
statevector(find(isExcluded(posdata(:,1), excludeperiods))) = -1; % Based on exclude time,

% spikes.data fields: time x y dir not_used amplitude(highest variance channel) posindex x-sm y-sm dir-sm
if size(spikes,2)>7
    posx = 8; posy = 9;
else
    posx = 2; posy = 3;
end
posx = 2; posy = 3;
posindexfield = 7;

if ~isempty(spikes)
    spikes = spikes(:,[1 posx posy posindexfield]);
    spikes(:,5) = statevector(spikes(:,4)); %add the traj state for each spike
    %posindexfield = isdatafield(spikesfields,'posindex');
    
else
    spikes = [0 0 -1];
end

trajnum = max(statevector);
timestep = statematrix.time(2,1) - statematrix.time(1,1);
goodspikes = [];
goodspikeind = (spikes(:,5) ~= -1);
%create a list of the non sharp-wave spikes [time traj linearloc]
goodspikes = spikes(goodspikeind,:);
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

%DR create seq indices for plotting w binsize>1
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
        
        %DR use edges for the whole track
        %         binx = output.allseqX1(:,2)';
        %         biny = output.allseqY1(:,2)';
        
        % Shantanu - I am using hist2 in ~/Src/Matlab/sj_Utility
        [output.occupancy{i} output.xticks{i} output.yticks{i}] = hist2(tmpposition(:,1), tmpposition(:,2), binx, biny);
        %DR sometimes hist2 messed up and addes an extra non-integer edge at the end of x edges... recrop matrix according to binx, biny
        output.occupancy{i} = output.occupancy{i}( 1:length(biny), 1:length(binx));
        output.xticks{i}  = output.xticks{i}(1,1:length(binx));
        output.yticks{i} = output.yticks{i}(1,1:length(biny));
        
        %DR create list as long as ticks for anchoring/plotting with a bin size more than 1
        output.seqX{i} = [output.allseqX(find(output.xticks{i}(1) == output.allseqX(:,2)),1):length(output.xticks{i})];
        output.seqY{i} = [output.allseqY(find(output.yticks{i}(1) == output.allseqY(:,2)),1):length(output.yticks{i})];
        
        if ~isempty(tmpspikes)
            
            %1) Get spikerate from spikes and occupancy
            [output.spikes{i} BX BY] = hist2(tmpspikes(:,1), tmpspikes(:,2), binx, biny);
            %DR sometimes hist2 messed up and addes an extra non-integer edge at the end of edges... recrop matrix according to binx, biny
            output.spikes{i} = output.spikes{i}(1:length(biny), 1:length(binx));
            BX  = BX(1,1:length(binx));
            BY = BY(1,1:length(biny));
            
            nonzero = find(output.occupancy{i} ~= 0);
            output.spikerate{i} = zeros(size(output.spikes{i}));
            skipp = 0;
            %             try
            output.spikerate{i}(nonzero) = output.spikes{i}(nonzero) ./(timestep* output.occupancy{i}(nonzero));
            %             catch
            %                 skipp = 1;
            %             end
            %             if skipp  == 1;
            %                 break
            %             end
            
            %2) Smooth spikerate and occupancy
            g = gaussian2(std,(6*std));
            output.smoothedspikerate{i} = filter2(g,(output.spikerate{i}));
            %smoothedoccupancy = [];
            %smoothedoccupancy = zeros(size(output(i).spikes));
            output.smoothedoccupancy{i} = filter2(g, output.occupancy{i});
            
            % 3) Turn occupancy to seconds and set spikerate wherever occupancy
            %is < threshold occupancy in seconds to 0
            
            output.occupancy{i} = timestep*output.occupancy{i};
            output.smoothedoccupancy{i} = timestep*output.smoothedoccupancy{i};
            
            %zero = find(smoothedoccupancy == 0);
            zero = find(output.smoothedoccupancy{i} <= threshocc);
            %output.smoothedspikerate{i}(zero) = -2; %no occupancy is negative and diff/darker from occupancy but no spikes
            try
                output.smoothedspikerate{i}(zero) = -1;
            catch
                output.smoothedspikerate{i} = []
            end
            
        else
            output.spikes{i} = [];
            output.spikerate{i} = [];
            output.smoothedspikerate{i} = [];
            output.smoothedoccupancy{i} = [];
        end
        
    else
        
        output.occupancy{i} = [];
        output.xticks{i}  = [];
        output.yticks{i} = [];
        output.spikes{i} = [];
        output.spikerate{i} = [];
        output.smoothedspikerate{i} = [];
        output.smoothedoccupancy{i} = [];
        
    end
    
    
    
end


%DR append a 2d map with all trajs
i = 5;
tmpposition = goodlocations((goodlocations(:,4) > 0 & goodlocations(:,4) < 5), [2 3]);
tmpspikes = goodspikes((goodspikes(:,5) > 0 & goodspikes(:,5) < 5),[2 3]);

% 
% tmpposition = goodlocations((goodlocations(:,4) == 1 | goodlocations(:,4) == 2 | goodlocations(:,4) == 3 | goodlocations(:,4) = 4) ,[2 3]);
% tmpspikes = goodspikes((goodspikes(:,5) == 1 || goodlocations(:,5) == 2 || goodlocations(:,5) == 3 || goodlocations(:,5) == 4),[2 3]);
if ~isempty(tmpposition)
    minx = floor(min(tmpposition(:,1)));   
    maxx = ceil(max(tmpposition(:,1)));
    binx = (minx:binsize:maxx);
    miny = floor(min(tmpposition(:,2)));
    maxy = ceil(max(tmpposition(:,2)));
    biny = (miny:binsize:maxy);
%     
%      [output.occupancy{i} output.xticks{i} output.yticks{i}] = hist2(tmpposition(:,1), tmpposition(:,2), binx, biny);  %this function is slow and buggy... just do 2d hist by hand:

xr = interp1(binx, 1:numel(binx), tmpposition(:,1), 'nearest');
yr = interp1(biny, 1:numel(biny), tmpposition(:,2), 'nearest');

xryr = [xr yr];
xryr = xryr(~any(isnan(xryr),2),:);

output.occupancy{i} = accumarray( [xryr] , 1, [numel(binx) numel(biny)]);


    
    
    
    
    %DR sometimes hist2 messed up and addes an extra non-integer edge at the end of x edges... recrop matrix according to binx, biny
    output.occupancy{i} = output.occupancy{i}( 1:length(biny), 1:length(binx));
    output.xticks{i}  = output.xticks{i}(1,1:length(binx));
    output.yticks{i} = output.yticks{i}(1,1:length(biny));
    output.seqX{i} = [output.allseqX(find(output.xticks{i}(1) == output.allseqX(:,2)),1):length(output.xticks{i})];
    output.seqY{i} = [output.allseqY(find(output.yticks{i}(1) == output.allseqY(:,2)),1):length(output.yticks{i})];
    if ~isempty(tmpspikes)
        
        %1) Get spikerate from spikes and occupancy
        [output.spikes{i} BX BY] = hist2(tmpspikes(:,1), tmpspikes(:,2), binx, biny);
        %DR sometimes hist2 messed up and addes an extra non-integer edge at the end of edges... recrop matrix according to binx, biny
        output.spikes{i} = output.spikes{i}(1:length(biny), 1:length(binx));
        BX  = BX(1,1:length(binx));
        BY = BY(1,1:length(biny));
        
        nonzero = find(output.occupancy{i} ~= 0);
        output.spikerate{i} = zeros(size(output.spikes{i}));
        skipp = 0;
        %             try
        output.spikerate{i}(nonzero) = output.spikes{i}(nonzero) ./(timestep* output.occupancy{i}(nonzero));
        %             catch
        %                 skipp = 1;
        %             end
        %             if skipp  == 1;
        %                 break
        %             end
        
        %2) Smooth spikerate and occupancy
        g = gaussian2(std,(6*std));
        output.smoothedspikerate{i} = filter2(g,(output.spikerate{i}));
        %smoothedoccupancy = [];
        %smoothedoccupancy = zeros(size(output(i).spikes));
        output.smoothedoccupancy{i} = filter2(g, output.occupancy{i});
        
        % 3) Turn occupancy to seconds and set spikerate wherever occupancy
        %is < threshold occupancy in seconds to 0
        
        output.occupancy{i} = timestep*output.occupancy{i};
        output.smoothedoccupancy{i} = timestep*output.smoothedoccupancy{i};
        
        %zero = find(smoothedoccupancy == 0);
        zero = find(output.smoothedoccupancy{i} <= threshocc);
        %output.smoothedspikerate{i}(zero) = -2; %no occupancy is negative and diff/darker from occupancy but no spikes
        try
            output.smoothedspikerate{i}(zero) = -1;
        catch
            output.smoothedspikerate{i} = []
        end
        
    else
        output.spikes{i} = [];
        output.spikerate{i} = [];
        output.smoothedspikerate{i} = [];
        output.smoothedoccupancy{i} = [];
    end
    
else
    
    output.occupancy{i} = [];
    output.xticks{i}  = [];
    output.yticks{i} = [];
    output.spikes{i} = [];
    output.spikerate{i} = [];
    output.smoothedspikerate{i} = [];
    output.smoothedoccupancy{i} = [];
end


if appendindex == 1
    output.index{1} = index;
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


