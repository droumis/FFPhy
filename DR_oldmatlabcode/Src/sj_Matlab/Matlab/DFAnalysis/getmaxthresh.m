function out = getmaxthresh(index, excludetimes, event, cellinfo, varargin)
% GETMAXTHRESH looks through each valid time segment and determines the
% maxthresh of events that occur during those time segments.
%
% index: [day epoch]
% event: event structure, such as ripple, gammal, or gammah

% Options:
% tetlist: A string that defines which tetrodes to use. Default is to
%   return maxthresh for all tetrodes with more than 1 cell.


tetlist = '$numcells > 1';

for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'tetlist'
                tetlist = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% assign a temporary variable for rip and define which tetrodes to use
e = event{index(1)}{index(2)};
c = cellfetch(cellinfo{index(1)}{index(2)},tetlist);
tetlist = unique(c.index(:,1));
clear event

% match max thresh values to included times
includetimes = getincludedtimes(excludetimes);
maxthresh = nan(size(includetimes,1),length(tetlist));
baseline = nan(length(tetlist),1);
stdeviation = nan(length(tetlist),1);
totalbaseline = nan(length(tetlist),1);
totalstd = nan(length(tetlist),1);

for t = 1:length(tetlist)
    try
        event = e{(tetlist(t))};
        idx = lookup(event.starttime, includetimes(:,1));
        maxthresh(idx,t) = event.maxthresh;
        baseline(t) = event.baseline;
        stdeviation(t) = event.std;
        totalbaseline(t) = event.totalbaseline;
        totalstd(t) = event.totalstd;
    end
end

out.maxthresh = maxthresh;
out.std = stdeviation;
out.baseline = baseline;
out.totalbaseline = totalbaseline;
out.totalstd = totalstd;
out.index = index;

end
