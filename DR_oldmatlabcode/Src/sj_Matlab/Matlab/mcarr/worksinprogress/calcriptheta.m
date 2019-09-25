function out = calcriptheta(index,excludetimes,eeg,ripples,cellinfo, varargin)

%  Computes the average CA1 and CA3 theta power and CA1-CA3 theta coherence
%  for the 400 ms before SWR detect and the correlation between theta power
%  and ripple power. Takes pairs of CA3 -CA1 tetrodes as input
%
%

%  Options:
%       fpass-- Determines the frequency range for computing coherence.
%           Default: [2 20]
%       average_trials-- Determines if events are averaged or not.
%           Default: 0
%       coherence_window-- Determines the sliding window used to compute
%           the event triggered coherence. Default: [0.4 0.4]
%       event_window--Determines the size of the window around each
%           triggering event. Default: [0.8 0.8]
%       cellfilter--Determines which tetrodes are used to detect triggering
%           events. Default: 'isequal($area,''CA1'')'

%parse the options
params = {};
params.Fs = 1500;
params.fpass = [2 20];
params.trialave = 0;
win = [0.8 0.8];
nonreference = 1;
windows = [0.4 0.4];
minthresh = 3;

for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})       
        switch(varargin{option})
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'aveargetrials'
                params.trialave = varargin{option+1};
            case 'event_window'
                win = varargin{option+1};
            case 'coherence_window'
                windows = varargin{option+1};
            case 'nonreference'
                nonreference = varargin{option+1};
            case 'minthresh'
                minthresh = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

out = [];

%Identify ripples, apply excludetimes
riptimes = getripples([index(1,1) index(1,2)], ripples, cellinfo, ...
    'cellfilter', '(isequal($area, ''CA1''))', 'excludeperiods', ...
    excludetimes,'minstd',minthresh);

if ~isempty(riptimes)
    
    %Exclude ripples occuring close together for preceding modulationg
    valid_ripples = [1000; riptimes(2:end,2)-riptimes(1:end-1,1)];
    valid_ripples = valid_ripples > 1;

    %Initialize eeg
    e1 = eeg{index(1,1)}{index(1,2)}{index(1,3)};
    e1times = geteegtimes(e1);

    e2 = eeg{index(1,1)}{index(1,2)}{index(1,4)};
    e2times = geteegtimes(e2);

    if length(e1times)>length(e2times)
        temp = lookup(e2times,e1times);
        e1times = e1times(temp);
        e1 = e1.data(temp);
        e2 = e2.data;
    elseif length(e2times)>length(e1times)
        temp = lookup(e1times,e2times);
        e1 = e1.data;
        e2 = e2.data(temp);
    elseif length(e1times)==length(e2times)
        e1 = e1.data;
        e2 = e2.data;
    end
    starttime = e1times(1);
    endtime = (length(e1)-1) * (1 / params.Fs);

    % Define triggering events as the start of each ripple
    triggers = riptimes(valid_ripples,1)-starttime;
    
    %Remove triggering events that are too close to the beginning or end
    while triggers(1)<win(1) + 0.5
        triggers(1) = [];
    end
    while triggers(end)> endtime-win(2)-0.5
        triggers(end) = [];
    end

    % Calculate the event triggered coherence
    data1 = createdatamatc(e1,triggers,params.Fs,[win(1) win(2)]);
    data2 = createdatamatc(e2,triggers,params.Fs,[win(1) win(2)]);

    [C,phi,S12,S1,S2,t,f] = cohgramc(data1,data2,[windows(1) windows(2)],params);

    %Determine theta coherence for 400ms before SWR
    time = lookup([-0.4 0],t-win(1));
    theta = lookup([6 12],f);
    
    out.coherence = squeeze(mean(mean(C(time(1):time(end),theta(1):theta(2),:),2),1));
    
    % Compute a z-scored spectrogram using the mean and std for the entire session
    P = mtspecgramc(e1,[windows(1) windows(2)],params);
    meanP = mean(P)';
    stdP = std(P)';
    
    out.ca1_power = squeeze(mean(S1(time(1):time(end),theta(1):theta(2),:),1));
    for i = 1:size(out.ca1_power,2)
        out.ca1_power(:,i) = (out.ca1_power(:,i)-meanP(theta(1):theta(2)))./stdP(theta(1):theta(2));
    end
    out.ca1_power = squeeze(mean(out.ca1_power,1));
    
    P = mtspecgramc(e2,[windows(1) windows(2)],params);
    meanP = mean(P)';
    stdP = std(P)';
    
    out.ca3_power = squeeze(mean(S2(time(1):time(end),theta(1):theta(2),:),1));
    for i = 1:size(out.ca3_power,2)
        out.ca3_power(:,i) = (out.ca3_power(:,i)-meanP(theta(1):theta(2)))./stdP(theta(1):theta(2));
    end
    out.ca3_power = squeeze(mean(out.ca3_power,1));
    
    out.starttime = triggers + starttime;
end
