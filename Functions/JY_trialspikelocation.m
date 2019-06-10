function [out] = JY_trialspikelocation(index, excludetimes, spikes,data,linpos, varargin)

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

spikesfields = spikes{index(1)}{index(2)}{index(3)}{index(4)}.fields;
spikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data;
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

if ~isempty(spikes)
    spikes = spikes(:,[1 posx posy posindexfield]); %columns: time, x y, posindex

else
    spikes = [0 0 -1];
end

output={};

% filter out excluded times

timestep = mean(diff(pos(:,1)));;
%filter out exclude times
indgoodspikes = ~isExcluded(spikes(:,1), excludetimes);
goodspikes = spikes(indgoodspikes,:);  %select spikes not excluded by exclude times
indgoodpos =  ~isExcluded(pos(:,1), excludetimes);
%goodpos = pos(indgoodpos,:);  %select spikes not excluded by exclude times

goodpos=pos;


% match spikes to each trial
trajectoryfiringrate=[];
trajectory={};
for triali=1:size(trajectorytimes,1)
    trialspikes=isExcluded(goodspikes(:,1)*10000, trajectorytimes(triali,:));
    trialposindex=isExcluded(goodpos(:,1)*10000, trajectorytimes(triali,:));
    trajectory{triali}.spikes=goodspikes(find(trialspikes==1),:);
    trajectory{triali}.pos=goodpos(trialposindex==1,:);
end
    
 
% match spikes to inter trial
intertrialtimes=[];
for ii=1:size(trajectorytimes,1)-1
    tstart=trajectorytimes(ii,2);
    tend=trajectorytimes(ii+1,1);
    intertrialtimes=[intertrialtimes; tstart tend];
end
    
% match spikes to each inter trial
intertrialfiringrate=[];
intertrial={};
for triali=1:size(intertrialtimes,1)
    trialspikes=isExcluded(goodspikes(:,1)*10000, intertrialtimes(triali,:));
    trialposindex=isExcluded(pos(:,1)*10000, intertrialtimes(triali,:));
    intertrial{triali}.spikes=goodspikes(find(trialspikes==1),:);
    intertrial{triali}.pos=goodpos(trialposindex==1,:);
end


out.trajectory=trajectory;
out.trajectorybarrier=trajectorybarrier;
out.intertrial=intertrial;
out.allpos=pos(:,2:3);
out.allspiketimes=goodspikes(:,1)*10000;

warning('ON','MATLAB:divideByZero');