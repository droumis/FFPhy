function [out] = getspikesperripple(index, excludeperiods, ripples, spikes, varargin)
% [out] = getspikesperripple(index, excludeperiods, ripples, spikes, varargin)
%
%   This function calculates the number of spikes each cell fired during
%   each ripple.
%
%   index [day epoch tetrode cell]
%
%   options are
%	'minthresh',
%		     specifies the minimum threshold of a valid ripple event
%   'appendindex' , 1 or 0, default 0
%           set to 1 to append the cell index to the output [day epoch
%           value]
%   'binsize', specifies the binsize for getincludedbinedges. If not
%           specified, binedges is not computed.
%   out = out.out   An R x C sized matrix where each entry is the number of
%                   spikes cell C fired during ripple R
%         out.index [D E T C], gives the identity of the cells for each
%                   column in out.out
%         out.times [starttime endtime], givest the starttime and endtime
%                   of the ripples for each row in out.out

% assign the options
minthresh = 0;
binsize = 0;
minduration = -1;

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'minthresh'
            minthresh = varargin{option+1};
        case 'binsize'
            binsize = varargin{option+1};
        case 'minduration'
            minduration = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%Use tetrode with most cells to detect ripples
maxcells = max(index(:,4));
ripindex = index(find(index(:,4)==maxcells),:);
rip = ripples{ripindex(1,1)}{ripindex(1,2)}{ripindex(1,3)};

%apply excludetimes
rvalidstart = ~isExcluded(rip.starttime, excludeperiods);
rvalidend = ~isExcluded(rip.endtime,excludeperiods);
rvalid = rvalidstart & rvalidend;

if (minthresh == 0)
    % get all the goodtimes
    rind = [rip.starttime(rvalid) rip.endtime(rvalid)];
else
    % get the indeces for the ripples with energy above minthreshold
    rvalidthresh = (rip.maxthresh > minthresh);
    rvalid = rvalid & rvalidthresh;
    rind = [rip.starttime(rvalid) rip.endtime(rvalid)];
end

%Count how many spikes occur during each ripple.
output = zeros(length(rind),length(index));

for i = 1:length(index)
    cell = index(i,:);
    spiketimes = spikes{cell(1,1)}{cell(1,2)}{cell(1,3)}{cell(1,4)}.data(:,1);
    %Step through each ripple
    for j = 1:length(rind)
        try
            output(j,i) = sum(spiketimes > rind(j,1) & spiketimes < rind(j,2) );
        catch
            output(j,i) = 0;
        end
    end
    
end

out.out = output;
out.index = index;
out.times = rind;

if binsize > 0
    if (minduration == -1)
        %Run getincludedbinedges
        binedges = getincludedbinedges(excludeperiods,binsize);
    else
        %Run getincludedbinedges with a defined minimum duration
        binedges = getincludedbinedges(excludeperiods,binsize,'minduration',minduration);
    end
    out.binedges = binedges;    
end

end