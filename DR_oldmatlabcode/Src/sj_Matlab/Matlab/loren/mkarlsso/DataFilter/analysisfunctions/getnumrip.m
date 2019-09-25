function [out] = getnumrip(index, excludeperiods, ripples, varargin)
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
%   'minstd', S
%           specifies the minimum ripple threshold (in standand deviations)
%           to count as a ripple
%   out is [rate proportiontime]
%       proprotiontime is proportion of included time during which ripples
%   were recorded
%       rate is number ripples/sec during included time

% assign the options
numtetrodes = 1;
minenergy = 0;
proptetrodes = [];
appendindex = 0;
minstd = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'numtetrodes'
            numtetrodes = varargin{option+1};
        case 'proptetrodes'
            proptetrodes = varargin{option+1};
        case 'appendindex'
            appendindex = varargin{option+1};
        case 'maxcell'
            maxcell = varargin{option+1};
        case 'minenergy'
            minenergy = varargin{option+1};
        case 'minstd'
            minstd = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

if ~isempty(proptetrodes)
    numtetrodes = round(size(index,1) * proptetrodes); %size(index,1) is total number of tetrodes input
end


%find tetrode with max number of cells per day
if maxcell == 1
    a = sortrows(epochs,1);
    days = unique(a(:,1));
    tetlistdays = [];
    for i = 1:length(days)  %find the tetrode with most cells on each day
        numcells = [0 0];
        d = days(i);
        eps = epochs(epochs(:,1)==d,2);
        for j = 1:length(eps)
            e = eps(j);
            if (~isempty(cellfilter))
                tetlist =  evaluatefilter(cellinfo{d}{e}, cellfilter);
                % get rid of the cell indeces and extract only the tetrode numbers
                tetlist = unique(tetlist(:,1))';
            end
            for s = 1:length(tetlist)
                tet = tetlist(s); %tetrode
                if ~exist('cellinfo')
                    cellinfo = loaddatastruct(animaldir, animalprefix, 'cellinfo');
                end
                if ismember(tet,numcells(:,1)) %if tetrode already listed in numcells
                    [row junk] = find(numcells(:,1) == tet);
                    numcells(row, 2) = numcells(row, 2) + length(cellinfo{d}{e}{tet});
                else
                    numcells = [numcells; tet length(cellinfo{d}{e}{tet}) ];
                end
            end
        end
        [y i] = max(numcells(:,2));
        tetlistdays = [tetlistdays; d numcells(i,1)];
    end
end

times = r.timerange(1):0.001:r.timerange(end);
nrip = zeros(size(times));
for t = 1:length(tetlist)
    tmprip = ripples{index(1,1)}{index(1,2)}{tetlist(t)};
    % get the indeces for the ripples with energy above minenergy
    rvalid = find((tmprip.energy > minenergy) & (tmprip.maxthresh >= minstd));
    rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
   
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
rip.nripples = nrip; %number of ripples on each tetrode

%apply excludetimes to nripples 
includetimes = ~isExcluded(rip.times, excludeperiods); %list of ones and zeros sampled every millisecond
if size(rip.nripples) ~= size(includetimes)
    includetimes = includetimes';
end
includerips = rip.nripples .* includetimes;

%calculate number of ripples per time 
a = zeros(1, length(includerips));
a(find(includerips >= numtetrodes)) = 1;
ripcount = length(find(diff(a) == 1));


%calculate ripple rate
if appendindex == 0
    out = ripcount;
elseif appendindex ==1
    out = [index ripcount];
end

end

