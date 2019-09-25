function out = plotthetamod(sind, tind, excludetimes, spikes, theta, varargin)
% out = plotthetamod(spike_index, theta_index, excludeperiods, spikes, theta, options)
% plots the theta moduation histogram for the specified cell based on the
% specified index into the theta structure.
% Excluded time periods are not included in calculation.
%
% Options:
%   'nbins', # of bins (default 12)
%   'appendind', 1 or 0 -- set to 1 to append the cell ind to the
%   output [tetrode cell value].  Default 0.
%

nbins = 12;
appendind = 0;
color = 'k';
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'color'
                color = varargin{option+1};
            case 'appendind'
                appendind = varargin{option+1};
            case 'nbins'
                nbins = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end


% get the spike times
s = spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.data(:,1);   %spike times

% get the eeg sample times
t = geteegtimes(theta{tind(1)}{tind(2)}{tind(3)});          % eeg clock time

tph = theta{tind(1)}{tind(2)}{tind(3)}.data(:,2);           % theta phase vector
sph = tph(lookup(s, t));                                    % spike phases



if ~isempty(excludetimes)
    totalexclude = sum(excludetimes(:,2) - excludetimes(:,1));
else
    totalexclude = 0;
end

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

sph = double(sph(goodspikes)) / 10000;

bins = -pi:(2*pi/nbins):pi;
count = histc(sph, bins);
out = bar(bins, count, 'hist')
set(out,'facecolor',color)
[m ph] = modulation(sph);
title(sprintf('cell %d %d %d %d, eegtet %d theta, mod %f', sind, tind(3), m));
pause
