function [out] = JY_alignspiketoreward(index, excludetimes, spikes,data, varargin)

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
reward = data{index(1)}{index(2)}.Events.Welltriggers;

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


%filter out exclude times
indgoodspikes = ~isExcluded(spikes(:,1), excludetimes);
goodspikes = spikes(indgoodspikes,:);  %select spikes not excluded by exclude times

% separate events for each well
spikeu=goodspikes(:,1)*10000;
uniquewells=unique(reward(:,2));
rewardwell={};
histogramedspikes={};
spiketimes={};

% time window
pret=4000; postt=6000;
binsize_plot=binsize; %% If you put binsize_plot=1000, then units are Nspikes/binsize, not inst. firing rate in Hz
timeaxis = -pret:binsize:postt;
bins_center = find(abs(timeaxis)<=100);
bins_resp = find(timeaxis>=0 & timeaxis<=100);
smwin=100; %Smoothing Window
binsize = 100;
nstd = round(binsize*2/binsize); g1 = gaussian(nstd, 1*nstd+1);


for i=1:size(uniquewells,1)
    rewardwell{i}=reward(reward(:,2)==uniquewells(i),:);

    cumspks=[];
  for ii=1:length(rewardwell{i})
        
        currtime = (rewardwell{i}(ii));
        currspks = spikeu(find( (spikeu>=(currtime-pret)) & (spikeu<=(currtime+postt)) ));
        %currspks = currspks-(currtime-pret); % Set to 0 at -pret
        %histspks = histc(currspks,[0:binsize:pret+postt]); % Make histogram
        currspks = currspks-(currtime); % Set to 0 at reward
        % Save spktimes
        cumspks=[cumspks; currspks];
        spiketimes{i}=cumspks;
        % Save histogram
        %histogramedspikes{i}{ii}=;
  end
        histspks = histc(spiketimes{i},[-pret:binsize:postt]); % Make histogram
        histspks = smoothvect(histspks, g1); % Smmoth histogram%histogramedspikes{i}{ii}=;
        histspkwell{i}=histspks;
end

out.spiketimes=spiketimes;
out.histogramespikes=histspkwell;
out.axis=-pret:binsize:postt;

warning('ON','MATLAB:divideByZero');