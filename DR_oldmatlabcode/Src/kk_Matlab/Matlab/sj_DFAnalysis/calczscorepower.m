function out = calczscorepower(index, excludetimes, eeg, varargin)
% function out = calczscorepower(index, excludetimes, eeg, varargin)
%
%  Returns the z-score power for the frequencies in fpass for all valid times.
%
%   out is a structure with the following fields:
%       zPower: This is the z-score power time series
%       time: This is the time series corresponding to zPower
%
%   Options:
%   fpass: Defines the frequency band of interest. Default is 65-90Hz.
%   win: Defines the window for which the spectrogram is calculated.

params = {};
params.Fs = 1500;

% set options
params.fpass = [65 90];
win = [0.5 0.5];

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'window'
                win = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% assign a temporary variable for eeg
e1 = eeg{index(1)}{index(2)}{index(3)};
etimes = geteegtimes(e1);
clear eeg

%Calculate the power for the band of interest across all time
[S,t,f] = mtspecgramc(e1.data,[win(1) win(2)],params);
S = 10*log10(S);
avgS = mean(S,2);

%apply excludetimes
t = t + e1.starttime;
includetimes = ~isExcluded(t, excludetimes);
avgS = avgS(includetimes);
t = t(includetimes);

%Compute the normalized time series
zS = (avgS -mean(avgS))./std(avgS);

out.zPower = zS;
out.time = t;

end