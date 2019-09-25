function [out] = DFAsj_getriprate(index, excludeperiods, ripples, varargin)
% out = getriptimes(index, excludeperiods, ripples, options)
%
%   index [day epoch]
%   tetlist is a list of tetrodes to include in the analysis
%
%   options are
%	'minenergy', E
%		     specifies the minimum energy of a valid ripple event
%   'minthresh', minthresh
%		     specifies a minimum threshold in stdev units for a valid
%			ripple event  (default 0)
%   'numtetrodes'
%           specifies number of tetrodes a ripple must be recorded on to be
%           included in analysis, default 1
%   'proptetrodes', examples: 1, 0.5, 0.25,
%           proportion of tetrodes a ripple must be recorded on to be
%           included in analysis
%   'appendindex' , 1 or 0, default 0
%           set to 1 to append the cell index to the output [day epoch
%           value]
%
%   out is [rate proportiontime]
%       proprotiontime is proportion of included time during which ripples
%   were recorded
%       rate is number ripples/sec during included time

% assign the options
numtetrodes = 1;
minenergy = 0;
minthresh = 0; % sd deviations
proptetrodes = [];
appendindex = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'numtetrodes'
            numtetrodes = varargin{option+1};
        case 'proptetrodes'
            proptetrodes = varargin{option+1};
        case 'appendindex'
            appendindex = varargin{option+1};
        case 'minenergy'
            minenergy = varargin{option+1};
        case 'minthresh'
            minthresh = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

if ~isempty(proptetrodes)
    numtetrodes = round(size(index,1) * proptetrodes); %size(index,1) is total number of tetrodes input
end

day = unique(index(:,1));
epoch = unique(index(:,2));


tetlist = index(:,3);
r = ripples{index(1,1)}{index(1,2)}{tetlist(1)}; % For time range for current epoch

times = r.timerange(1):0.001:r.timerange(end);
nrip = zeros(size(times));
nrip_size = zeros(length(tetlist),length(times));

for t = 1:length(tetlist)
    tmprip = ripples{index(1,1)}{index(1,2)}{tetlist(t)};
    % Shantanu - change from minenergy to minthresh
    if (minthresh == 0)
        % get all the times
        rtimes = [tmprip.starttime tmprip.endtime];
        rvalid = 1:length(tmprip.maxthresh); % All ripples are valid
    else
        % get the indices for the ripples with thresh above minthresh
        rvalid = find(tmprip.maxthresh > minthresh);
        rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
    end
    % For ripsize
    currsize=tmprip.maxthresh(rvalid); % Size of each ripple on current tetrode
    
    if ~isempty(rtimes)
        % create another parallel vector with bordering times for zeros
        nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];
        rtimes = reshape(rtimes', length(rtimes(:)), 1);
        rtimes(:,2) = 1;
        
        % For ripsize
        % Make currsize match the format of rtimes and add it there
        currsize(:,2) = currsize;
        currsize = reshape(currsize', length(currsize(:)), 1);
        rtimes(:,3) = currsize;
        %     rtimes_size = rtimes(:,1);
        %     rtimes_size(:,2) = currsize;
        
        % Update nrtimes
        nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
            length(nrtimes(:)), 1) ; r.timerange(2)];
        nrtimes(:,2) = 0;
        nrtimes(:,3) = 0;
        
        % Create a new list with all of the times in it with 1 or 0 for ripple
        % or not
        tlist = sortrows([rtimes ; nrtimes]);
        % Alternate method
        %     % create a new list with all of the times in it with ripsize instead of
        %     % 1 or 0
        %     tlist_size = sortrows([rtimes_size ; nrtimes]);
        
        % use interp to create a set of ones and zeros for each time
        % and add to nrip to get a cumulative count of the number of
        % ripples per timestep
        try
            nrip = nrip + interp1(tlist(:,1), tlist(:,2), times, 'nearest');
        catch
            keyboard
        end
        
        try
            nrip_size(t,:) = interp1(tlist(:,1), tlist(:,3), times, 'nearest');
        catch
            keyboard
        end
        
    end % if ~isempty riptime
    
    % mean, std and threshold of ripple envelope used for each tet in current epoch
    % Return a mean acroos all tets for comparing acroos tets later
    % Also, sj_calcriprate can figure out which tet is assigning ripsize,
    % and return its baseline, std and threshold. % Not implented now
    rip.baseline(t) = tmprip.baseline;
    rip.std(t) = tmprip.std;
    rip.thresh(t) = tmprip.threshold;
end
rip.times = times; %sampled every millisecond
clear times;
rip.nripples = nrip; %number of ripples on all tetrode
rip.ripsize = nrip_size; %ripsize per ms on each tet in diff columns: same size is reported from start to end of ripple


%calculate ripple rate
% No passing of stimtimes to sj_calcriprate unlike DFAsj_getriprate_noDIO

if appendindex == 0
    out = sj_calcriprate(rip, excludeperiods, numtetrodes);
elseif appendindex ==1
    out = sj_calcriprate(rip, excludeperiods, numtetrodes);
    out.index = index;
end

% Also, return mean, std and threshold of ripple envelope used for
% current epoch
out.rip_baseline = mean(rip.baseline);
out.rip_std = mean(rip.std);
out.rip_thresh = mean(rip.thresh);
