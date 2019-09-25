function out = calceegmod(sind, eind, excludetimes, eeg, spikes, varargin)
% out = calceegmod(spike_index, theta_index, excludeperiods, spikes, theta, options)
% plots the eeg moduation histogram for the specified cell based on the
% specified index into the eeg structure.
% Excluded time periods are not included in calculation.
%
% Options:
% 'nbins', # of bins (default 12)
% 'appendind', 1 or 0 -- set to 1 to append the cell ind to the
%   output [tetrode cell value].  Default 0.
% 'frequency,' to bandpass filter the eeg input enter the low and high
%   cutoff. Default is to calculate modulation based on phase of eeg input
%   structure. Must use 'frequency' option if feeding in raw eeg structure.
%

frequency = '';
appendind = 0;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendind = varargin{option+1};
            case 'nbins'
                nbins = varargin{option+1};
            case 'frequency'
                frequency = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% Filter eeg if necessary and determine eeg phase for all times
if ~isempty(frequency)
    e = eeg{tind(1)}{eind(2)}{eind(3)}.data;
    phase = angle(hilbert(eegfilt(e',1500,frequency(1),frequency(2))));
else
    try
        phase = double(eeg{eind(1)}{eind(2)}{eind(3)}.data(:,2));
    catch
        warning('EEG structure is not in correct form')
    end
end

% get the spike and eeg times
s = spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.data(:,1);
t = geteegtimes(eeg{eind(1)}{eind(2)}{eind(3)});

% lookup the phase for each spike
sph = phase(lookup(s,t));

if ~isempty(s)
    if ~isempty(sph)
        goodspikes = ~isExcluded(s, excludetimes);
    else
        goodspikes = 0;
    end
else
    goodspikes = nan;
end

numgoodspikes = sum(goodspikes);

sph = sph(goodspikes)./10000;

out.phase = sph;
if appendind
    out.index = sind;
end
end