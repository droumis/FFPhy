function out = plotspectrogram(index, excludetimes, eeg, varargin)
% function out = plotspectrogram(index, excludetimes, eeg, varargin)
%
%  Plots the spectrogram for a given eeg tetrode. If you use a time filter,
%  excluded times are set to zero in the spectrogram (dark blue vertical
%  lines).
%
%   out is a structure with the following fields:
%       fullspectrum-- This is the spectrum for the eeg tetrode ignoring
%           your excludetimes
%        spectrum-- This is the spectrum with excluded times blocked out.
%       time-- Time vector
%       frequency-- Frequency vector
%       index-- Only if appendindex is set to 1

params = {};
params.Fs = 1500;

% set options
params.fpass = [0 150];
appendindex = 0;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'appendindex'
                appendindex = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% assign a temporary variable for eeg
e1 = eeg{index(1)}{index(2)}{index(3)};

clear eeg;

% compute zscored full spectrum
[S,t,f]= mtspecgramc(e1.data,[10 0.05],params);

%apply excludetimes
includetimes = ~isExcluded(t+e1.starttime, excludetimes);
temp = (includetimes==0);
tempspec = S;
tempspec = tempspec(temp,:);

out.fullspectrum = S;
out.spectrum = tempspec;
out.time = t+e1.starttime;
out.frequency = f;

if (appendindex)
    out.index = index;
end
