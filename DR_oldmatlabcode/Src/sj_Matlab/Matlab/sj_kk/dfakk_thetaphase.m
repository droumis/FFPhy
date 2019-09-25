function out = dfakk_thetaphase(sind, tind, excludetimes, spikes, theta, varargin)
% out = plotthetamod(spike_index, theta_index, excludeperiods, spikes, theta, options)
% plots the theta moduation histogram for the specified cell based on the
% specified index into the theta structure.
% Excluded time periods are not included in calculation.
%
% Options:
%   'isi', burst- and nonburst- associated classification of spikes
        % specify NEGATIVE isi for burst-associated, POSITIVE for non-burst
%   'nbins', # of bins (default 12)
%   'appendind', 1 or 0 -- set to 1 to append the cell ind to the
%   output [tetrode cell value].  Default 0.
%

nbins = 12;
appendind = 0;
isi=0;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'isi'      % see below         
                isi = varargin{option+1};
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

% if isi specified, then filter spikes by isi
    % select spikes that have < isi on either side (Mizuseki--Buzsaki-2009)

    if isi ~= 0
        
        filteredspikes = [];
        
        if isi < 0              % filtering for burst-associated spikes
            s = [0 ; s ; inf];
            for k=2:(length(s)-1)
                if (abs(s(k)-s(k-1)) < -isi/1000) || (abs(s(k)-s(k+1)) < -isi/1000)
                    filteredspikes = [ filteredspikes ; s(k)];
                end
            end
            s = filteredspikes;
        elseif isi > 0          % filtering for nonburst-associated spikes
            s = [0 ; s ; inf];
            for k=2:(length(s)-1)
                if (abs(s(k)-s(k-1)) > isi/1000) && (abs(s(k)-s(k+1)) > isi/1000)
                    filteredspikes = [ filteredspikes ; s(k)];
                end
            end
            s = filteredspikes;            
        end
    end


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

if ~isempty(sph)     
    sph = double(sph(goodspikes)) / 10000;
    out.sph = sph;
    out.index = sind;
else                   % burst filtering may eliminate ALL spikes in epoch
    out.sph = zeros(0,1);
    out.index = sind;
end

return





out.index = sind;
out.tetindex = tind;
out.sph = sph;
out.Nspikes = numgoodspikes;
% Output stats also, although you will have to redo this after combining epochs
% Rayleigh test
out.stats = stats;
out.modln = m;
out.phdeg = phdeg;
% Von Mises Fit
out.kappa = kappa;
out.thetahat_deg = thetahat_deg;
out.prayl = prayl;
out.zrayl = zrayl;











