function [out] = JY_calctrajfiringrate(index, excludetimes, spikes,data,linpos, varargin)
%[out] = openfieldoccupancy(index, excludetimes, spikes,pos, options)
%
%Calculates the 2d occupancy normalized firing rate for the cell.  Does NOT
%separate by trajectory (to do so, see twodoccupancy)
%
% excludetime -- times want to exclude from analysis
% spikes - the 'spikes' cell array for the day you are analyzing
% pos - the output of nspike_fixpos
% index - [day epoch tetrode cell]
% options:
%  binsize- the length of each spatial bin (default 1 cm)
%  std - defines the shape of the 2d gaussian used to smooth spikerate.
%              (default 1)
%  appendindex - 0 or 1, 1 appends index infront of output (default 0 )
%
%The output is a structure with fields: occupancy, bin vector x (.xticks), bin vector
%y (.yticks), bin spike count (.spikes), occ normailized firing per bin (.spikerate), and smoothed occ
% normalized firing (.smoothedspikerate).


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
velocity = data{index(1)}{index(2)}.Pos.correcteddata(:,5);
trajectorytimes=data{index(1)}{index(2)}.Run(:,3:4);
trajectoryduration=data{index(1)}{index(2)}.Run(:,5);
trajectorysegments=cellfun(@(x) size(x,1),linpos{index(1)}{index(2)}.trialsegments(2,:))';
trajectorybarrier=data{index(1)}{index(2)}.Run(:,6);

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
goodpos = pos(indgoodpos,:);  %select spikes not excluded by exclude times
goodvelocity = velocity(unique(goodspikes(:,4)));
tmpposition = (goodpos(:,[2 3])); %xy position
tmpspikes = (goodspikes(:,[2 3]));


% match spikes to each trajectory
trajectoryfiringrate=[];
trajectoryspikes={};
for triali=1:size(trajectorytimes,1)
    trialspikes=isExcluded(goodspikes(:,1)*10000, trajectorytimes(triali,:));
    trajectoryfiringrate=[trajectoryfiringrate; sum(trialspikes)./(trajectoryduration(triali)/10000)];
    trajectoryspikes{triali}=goodspikes(trialspikes==1,:);
end
 

% calculate inter trial firing rate

% match spikes to inter trial
intertrialtimes=[];
for ii=1:size(trajectorytimes,1)-1
    tstart=trajectorytimes(ii,2);
    tend=trajectorytimes(ii+1,1);
    intertrialbarrier=trajectorybarrier(ii+1,1);
    intertrialtimes=[intertrialtimes; tstart tend intertrialbarrier];
end
    
% match spikes to each inter trial

intertrialfiringrate=[];
intertrial={};
for triali=1:size(intertrialtimes,1)
    trialspikes=isExcluded(goodspikes(:,1)*10000, intertrialtimes(triali,:));
    trialposindex=isExcluded(pos(:,1)*10000, intertrialtimes(triali,:));
    intertrialfiringrate=[intertrialfiringrate; (sum(trialspikes)./((intertrialtimes(triali,2)-intertrialtimes(triali,1))))/10000];
    intertrialspikes{triali}=goodspikes(trialspikes==1,:);
end


        
    
out.trajectoryfiringrate=trajectoryfiringrate;
out.intertrialfiringrate=intertrialfiringrate;

out.trajectoryspikes=trajectoryspikes;
out.intertrialspikes=intertrialspikes;

out.allspikes=goodspikes;
out.velocity=goodvelocity;
out.allvelocity=pos(:,5);
out.trajectorysegments=trajectorysegments;
out.trajectoryduration=trajectoryduration./10000;
out.trajectorybarrier=trajectorybarrier;
out.intertrialbarrier=intertrialtimes(:,3);

warning('ON','MATLAB:divideByZero');