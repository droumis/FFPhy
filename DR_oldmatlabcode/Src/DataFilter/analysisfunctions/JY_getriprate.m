function [out] = getriprate(index, excludeperiods, ripples, varargin)
% out = getriptimes(index, excludeperiods, ripples, options)
%
%   index [day epoch]
%   tetlist is a list of tetrodes to include in the analysis
%
%   options are
%	'minenergy', E
%		     specifies the minimum energy of a valid ripple event
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
minenergy = 5;
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
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

if ~isempty(proptetrodes)
    numtetrodes = round(size(index,1) * proptetrodes); %size(index,1) is total number of tetrodes input
end

tetlist = index(:,3);
r = ripples{index(1,1)}{index(1,2)}{tetlist(1)};

binsize=0.1;

times = r.timerange(1):binsize:r.timerange(end);
nrip = zeros(size(times));
for t = 1:length(tetlist)
    tmprip = ripples{index(1,1)}{index(1,2)}{tetlist(t)};
    if (minenergy == 0)
        % get all the times
        rtimes = [tmprip.starttime tmprip.endtime];
    else
        % get the indeces for the ripples with energy above minenergy
        rvalid = find(tmprip.energy > minenergy);
        rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
    end
    % create another parallel vector with bordering times for zeros
    nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];
    rtimes = reshape(rtimes', length(rtimes(:)), 1);
    rtimes(:,2) = 1;
    nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
        length(nrtimes(:)), 1) ; r.timerange(2)];
    nrtimes(:,2) = 0;
    % create a new list with all of the times in it
    tlist = sortrows([rtimes ; nrtimes]);
    % use interp to create a set of ones and zeros for each time
    % and add to nrip to get a cumulative count of the number of
    % ripples per timestep
    try
        nrip = nrip + interp1(tlist(:,1), tlist(:,2), times, 'nearest');
    catch
        keyboard
    end
end
rip.times = times;  %same sampling as pos, spikes, etc
clear times;
rip.nripples = nrip; %number of ripples active at a given times (tetrodes)

%calculate ripple rate
if appendindex == 0
    out.rate = calcriprate(rip, excludeperiods, numtetrodes);
    out.times=rip.times;
    out.nripples=rip.nripples;
elseif appendindex ==1
    tmpout = calcriprate(rip, excludeperiods, numtetrodes);
    out.rate = [index tmpout];
end

end
