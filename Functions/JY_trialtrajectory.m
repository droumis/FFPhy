function [out] = JY_trialtrajectory(index, excludetimes,data,linpos, varargin)

%
% returns the spikes for each trial and each intertrial interval (reward
% sites)
%
% excludetime -- times want to exclude from analysis
% spikes - the 'spikes' cell array for the day you are analyzing
% pos - the output of nspike_fixpos
% index - [day epoch tetrode cell]



appendindex = 0;
std = 1;
binsize = 1;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'std'
                std = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end
warning('OFF','MATLAB:divideByZero');

% get the data


pos = data{index(1)}{index(2)}.Pos.correcteddata;
trajectorybarrier=data{index(1)}{index(2)}.Run(:,6);
trajectorytimes=data{index(1)}{index(2)}.Run(:,3:4);

if (nargin < 6)
    std = 1;
else
    %std = user defined;
end

if (nargin < 5)
    binsize = 1;
else
    %binsize = user defined;
end

posx = 2;
posy = 3;
posindexfield = 7;


output={};

% filter out excluded times


goodpos=pos;


% trial trajectory

trajectory={};
for triali=1:size(trajectorytimes,1)
 
    trialposindex=isExcluded(goodpos(:,1)*10000, trajectorytimes(triali,:));
trajectory{triali}.pos=goodpos(trialposindex==1,:);
end
    




out.trajectory=trajectory;
out.trajectorybarrier=trajectorybarrier;

out.allpos=pos(:,2:3);


warning('ON','MATLAB:divideByZero');