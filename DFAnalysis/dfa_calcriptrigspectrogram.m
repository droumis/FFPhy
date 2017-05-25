function out = dfa_calcriptrigspectrogram(index, excludeperiods, eeg, events, varargin)
% out = calcripspectrum(index, excludeperiods, eeg,ripples,cellinfo, options)
%  Computes the spectrogram around the middle of each decoded event.
%  Options:
%       appendindex-- Determines whether index is included in output vector
%           Default: 1
%       fpass-- Determines the frequency range for computing spectrum.
%           Default: [2 350]
%       average_trials-- Determines if events are averaged or not.
%           Default: 0
%       spectrum_window-- Determines the sliding window used to compute
%           the event triggered spectrogram. Default: [0.1 0.01]
%       event_window--Determines the size of the window around each
%           triggering event. Default: [0.2 0.2]
%       cellfilter--Determines which tetrodes to use for ripple extraction.
%           Default is 'isequal($area, ''CA1'') & numcells>1)'
%  out is a structure with the following fields
%       S-- This is a MxNxT matrix of the spectrogram for each tetrode
%           M is time relative to triggering event, N is frequency, T is event
%       F-- Frequency vector
%       T-- time relative to triggering event
%       fit-- This is the fit based on the spectrum computed for the entire
%           epoch to normalize S. To reconstruct S without normalization,
%           add log10(frequency)*fit(2)+fit(1) to S.
%       index-- Only if appendindex is set to 1 (default)


%parse the options
params = {};
params.Fs = 1500;
params.fpass = [2 350];
params.trialave = 0;
params.tapers = [3 5];
win = [0.3 0.3];
cwin = [0.1 0.01];
appendindex = 1;
cellfilter = [];

for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'aveargetrials'
                params.trialave = varargin{option+1};
            case 'spectrum_window'
                cwin = varargin{option+1};
            case 'event_window'
                win = varargin{option+1};
            case 'cellfilter'
                cellfilter = varargin{option+1};
                 case 'eventtype'
                eventtype = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

  
if strcmp(eventtype,'ripples')
    eventtimes = events{index(1)}{index(2)}{index(3)}.starttime;
else % ripplecons, ripplekons, etc
    eventtimes = events{index(1)}{index(2)}{1}.starttime;
end
    eventtimes = eventtimes(~isExcluded(eventtimes,excludeperiods));


% Define EEG
e = eeg{index(1)}{index(2)}{index(3)}.data';
starttime = eeg{index(1)}{index(2)}{index(3)}.starttime;
endtime = (length(e)-1) * (1 / params.Fs);
clear eeg ripple cellinfo

% Define triggering events as the start of each ripple
starttime = double(starttime);
triggers = eventtimes(:,1)-starttime;

%Remove triggering events that are too close to the beginning or end
while triggers(1)<win(1)
    triggers(1) = [];
end
while triggers(end)> endtime-win(2)
    triggers(end) = [];
end

% Calculate the event triggered spectrogram
[S,t,f] = mtspecgramtrigc(e,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);
out.triggers = triggers;

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

out.S = S;
out.t = t - win(1);
out.f = f;

    
if appendindex
	out.index = index;
end
    
end