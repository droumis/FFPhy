function out = dfakk_sta(sind, eind, excludetimes, spikes, eeg, varargin)
% out = dfakk_sta(spike_index, eeg_index, excludeperiods, spikes, theta, options)

% outputs eeg windows surrounding spikes

% Excluded time periods are not included in calculation.
%
% Options:
%   'isi', burst- and nonburst- associated classification of spikes
        % specify NEGATIVE isi for burst-associated, POSITIVE for non-burst
%   'nbins', # of bins (default 12)

%   'appendind', 1 or 0 -- set to 1 to append the cell ind to the
%   output [tetrode cell value].  Default 0.

%   'window' in seconds -- size of eeg window around spike
            %  ex. [0.5 0.5]
%

appendind = 0;
isi = 0;
window = [0.5 0.5];

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'isi'      % see below         
                isi = varargin{option+1};
            case 'appendind'
                appendind = varargin{option+1};
            case 'window'
                window = varargin{option+1};
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

% if isi argument is specified, then filter spikes by isi
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
    
% remove spikes that fall in excludetimes    

if ~isempty(s)
   goodspikes = ~isExcluded(s, excludetimes);
   s = s(goodspikes);
end

% retrieve eeg sample times
eegdata = eeg{eind(1)}{eind(2)}{eind(3)};
t = geteegtimes(eegdata);          % eeg clock time
Fs = round(eegdata.samprate);
window_presamp = Fs*window(1);
window_postsamp = Fs*window(2);

%initialize output
out.eegwindows = [];
out.index = sind;

if ~isempty(s)     
    for k=1:length(s)
        centerindex = lookup(s(k),t);
        if (s(k)<(t(end)-window(2))) && (s(k)>(t(1)+window(1)))      % ignore spikes that occur at very beginning and end of epoch
            out.eegwindows = [out.eegwindows ; eegdata.data((centerindex-window_presamp):(centerindex+window_postsamp))'];
        end
    end
else                   % possibly no spikes in this epoch
    out.eegwindows = [];
end

return
















