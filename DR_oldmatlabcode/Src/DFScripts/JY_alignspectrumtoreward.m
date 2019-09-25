function [out] = JY_alignspectrumtoreward(index, excludetimes,data,eeg,meaneegspectrograms,tetinfo, varargin)

% Based on Maggie's calcriptriggeredspectrogram
% get the times of reward events from data
% excludetime -- times want to exclude from analysis
%

% index - [day epoch tetrode cell]



appendindex = 0;

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
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

reward = [data{index(1)}{index(2)}.Run(1:end-1,2) data{index(1)}{index(2)}.Run(1:end-1,4)];

uniquereward=unique(reward(:,1));

%for ii=1:size(uniquereward,1)
    reward1=reward(reward(:,1)==uniquereward(1),:);
%end
   
    

% separate events for each well

uniquewells=unique(reward(:,2));
rewardwell={};


for i=1:size(uniquewells,1)
    rewardwell{i}=reward(reward(:,2)==uniquewells(i),:);
end

out.rewardedevents=rewardwell;

reward=reward1;

%% generate spectrogram

%parse the options
params = {};
params.Fs = 1500;
params.fpass = [1 250];
params.trialave = 0;
cwin = [250 10]/1000;
win = [500 2000]/1000;
appendindex = 1;
cellfilter = [];
params.tapers = [3 5];



% Define EEG
e = meaneegspectrograms{index(1)}{index(2)}{index(4)}.unreferenced';
%e = eeg{index(1)}{index(2)}{index(4)}.data';
starttime = meaneegspectrograms{index(1)}{index(2)}{index(4)}.starttime;
%starttime = eeg{index(1)}{index(2)}{index(4)}.starttime;
endtime = (length(e)-1) * (1 / params.Fs);
clear eeg ripple cellinfo

% Define triggering events as the start of each ripple
triggers = reward(:,2)/10000-starttime;

%Remove triggering events that are too close to the beginning or end
while triggers(1)<win(1)
    triggers(1) = [];
end
while triggers(end)> endtime-win(2)
    triggers(end) = [];
end

% Calculate the event triggered spectrogram
[S,t,f] = mtspecgramtrigc(e,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);

% Compute a z-scored spectrogram using the mean and std for the entire session
P = mtspecgramc(e,[cwin(1) cwin(1)],params);
clear e triggers riptimes
meanP = mean(P);
stdP = std(P);
clear P
for i = 1:size(S,1)
    for j = 1:size(S,3)
        S(i,:,j) = (S(i,:,j) - meanP)./stdP;
    end
end

% for ii=1:5
%     
%     figure;
%     imagesc(t-win(1),f,S(:,:,ii)');
%     %imagesc(tmt-win(1),fmt,epocheegmean',[0, 20]);
%     set(gca,'YDir','normal');
%     string = sprintf('%d',t);
%     %title(string);
%     colormap('jet')
%     colorbar
%     close
% end

out.epochmeanzscore = mean(S,3);
out.times = t - win(1);
out.freq = f;

if appendindex
    out.index = index;
end
warning('ON','MATLAB:divideByZero');
end

