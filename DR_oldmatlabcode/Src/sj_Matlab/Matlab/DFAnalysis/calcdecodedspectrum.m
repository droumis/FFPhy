function out = calcdecodedspectrum(index, excludetimes, decodetimes, eeg,varargin)

% out = calcdecodespectrum(index, excludeperiods, decodetimes, eeg, options)
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
win = [0.2 0.2];
cwin = [0.1 0.01];
appendindex = 1;

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
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

if ~isempty(decodetimes.eventtime)
    % Define EEG
    e = eeg{index(1)}{index(2)}{index(3)}.data';
    starttime = eeg{index(1)}{index(2)}{index(3)}.starttime;


    % Define triggering events as the start of each decoded event
    triggers = decodetimes.eventtime(:,1)-starttime;

    % Calculate the event triggered spectrogram
    [S,t,f] = mtspecgramtrigc(e,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);

    % Normalize the spectrogram using the average spectrum for the session
    P = mtspecgramc(e,[cwin(1) cwin(2)],params);
    meanP = mean(P);
    stdP = std(P);
    for i = 1:size(S,1)
        for j = 1:size(S,3)
            S(i,:,j) = (S(i,:,j) - meanP)./stdP;
        end
    end

    out.S = S;
    out.t = t - win(1);
    out.f = f;
    out.mean = meanP;
    out.std = stdP;
    
    if appendindex
        out.index = index;
    end
    
else
    out.S = [];
    out.t = [];
    out.f = [];
    out.mean = [];
    out.std = [];
    if appendindex
        out.index = index;
    end
end

end