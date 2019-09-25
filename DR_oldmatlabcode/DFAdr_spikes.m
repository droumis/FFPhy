function out = DFAdr_spikes(index, excludeperiods, spikes, linpos, pos, varargin)

%output spike data and ISI

% parse the options
appendindex = 1;
binsize = 1; % cm square 
std = 2; 
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

% Spikes
spikestruct = spikes{index(1)}{index(2)}{index(3)}{index(4)};
spikes = spikestruct.data;

% Position Data
posdata = pos{index(1)}{index(2)}.data;

statematrix = linpos{index(1)}{index(2)}.statematrix;
statevector = statematrix.traj;

% Use Exclude Periods for TimeFilter version in addition to statevector=-1
% that already exists based on behavestate
% statematrix.time and posdata(:,1) are equivalent
statevector(find(isExcluded(posdata(:,1), excludeperiods))) = -1; % Based on exclude time, 

% spikes.data fields: time x y dir not_used amplitude(highest variance channel) posindex x-sm y-sm dir-sm
posx = 2; posy = 3;
posindexfield = 7;

if ~isempty(spikes)
    spikes = spikes(:,[1 posx posy posindexfield]);
    % Find any posidxs that are outside the range of linpos (statevector)
    outliers = find(spikes(:,4) > length(statevector));
    spikes(outliers,:)=[];
    spikes(:,5) = statevector(spikes(:,4)); %add the traj state for each spike
else
    spikes = [0 0 -1];
end


timestep = posdata(2,1) - posdata(1,1);

%index(3),index(4), size(spikes)

goodspikes = [];
goodspikes = spikes;
% Non sharp-wave spikes
goodspikeind = (spikes(:,5) ~= -1);
goodspikes = spikes(goodspikeind,:);

%output the goodspikes, along with info from the spikestruct
% output.index = spikestruct;
output.goodspikes = goodspikes;

if appendindex == 0
    out = output;
elseif appendindex == 1
    output.index = index;
    out = output;
end

warning('ON','MATLAB:divideByZero');


