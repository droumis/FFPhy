function out = dfakk_thetaphase(sind, tind, excludetimes, spikes, theta, varargin)

%  singlecelleeganal
%  retrieves theta phases of spikes
%   'nbins', # of bins (default 12)

nbins = 24;
isi = 0;

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'isi'      % see below         
                isi = varargin{option+1};
            case 'nbins'
                nbins = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% initialize output
out.index = sind;
out.sph = nan;
out.included_duration = nan;

% epoch duration + included_duration
starttime_epoch = spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.timerange(1)/10000;
endtime_epoch = spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.timerange(2)/10000;
    epoch_duration = endtime_epoch - starttime_epoch;
% calculate included duration
if ~isempty(excludetimes)
    exclude_duration = sum(excludetimes(:,2) - excludetimes(:,1));
else
    exclude_duration = 0;
end
included_duration = epoch_duration - exclude_duration;

% get the spike times
if ~isempty(spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.data)
    s = spikes{sind(1)}{sind(2)}{sind(3)}{sind(4)}.data(:,1);
    out.included_duration = included_duration;
else    % empty cluster
    out.sph = zeros(0,1);
    out.included_duration = included_duration;
    return
end

% filter for spikes in included periods
goodspikes = ~isExcluded(s, excludetimes);

if sum(goodspikes) > 0
    % get their theta phases
    thetatimevec = geteegtimes(theta{tind(1)}{tind(2)}{tind(3)});        
        inds = lookup(s(goodspikes), thetatimevec);
    tph = theta{tind(1)}{tind(2)}{tind(3)}.data(:,2);                   % theta phase vector
    sph = tph(inds);                                         % spike theta phases
    % output
    out.sph = double(sph)/10000;
else
    out.sph = zeros(0,1);    
end

    

return








% % if isi specified, then filter spikes by isi
%     % select spikes that have < isi on either side (Mizuseki--Buzsaki-2009)
%     
%     if isi ~= 0
%         
%         filteredspikes = [];
%         
%         if isi < 0              % filtering for burst-associated spikes
%             s = [0 ; s ; inf];
%             for k=2:(length(s)-1)
%                 if (abs(s(k)-s(k-1)) < -isi/1000) || (abs(s(k)-s(k+1)) < -isi/1000)
%                     filteredspikes = [ filteredspikes ; s(k)];
%                 end
%             end
%             s = filteredspikes;
%         elseif isi > 0          % filtering for nonburst-associated spikes
%             s = [0 ; s ; inf];
%             for k=2:(length(s)-1)
%                 if (abs(s(k)-s(k-1)) > isi/1000) && (abs(s(k)-s(k+1)) > isi/1000)
%                     filteredspikes = [ filteredspikes ; s(k)];
%                 end
%             end
%             s = filteredspikes;
%         end
%     end


