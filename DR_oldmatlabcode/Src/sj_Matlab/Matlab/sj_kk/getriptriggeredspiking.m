function [out] = getriptriggeredspiking(index, excludeperiods, ripples, spikes, varargin)
% [out] = getriptriggeredspiking(index, excludeperiods, ripples, spikes, varargin)
%
%   This function creates a ripple triggered spiking average. The ripples
%   are defined by the first group of cells and the spikes are from the
%   second group of cells.
%
%   index [day epoch tetrode cell tetrode cell]
%
%   options are
%	'minthresh',
%		     specifies the minimum threshold of a valid ripple event
%   'window', 1x2 vector specifies the window before and after each included ripple.
%                   Default is 100 mseconds before and 15 seconds after
%                   ripple start time.
%   out = out.out   An R x C sized matrix where each entry is the number of
%                   spikes cell C fired during ripple R
%         out.index [D E T C], gives the identity of the cells for each
%                   column in out.out
%         out.times [starttime endtime], givest the starttime and endtime
%                   of the ripples for each row in out.out

% assign the options
minthresh = 0;
window = [0.1 0.1];
binsize = 0.005;

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'minthresh'
            minthresh = varargin{option+1};
        case 'window'
            window = varargin{option+1};
        case 'binsize'
            binsize = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end


%Create times of all ripples using the first group of tetrodes

tetlist = index(:,3);                                       % tetlist: the "trigger" tetrodes
tetlist = unique(tetlist);
r = ripples{index(1,1)}{index(1,2)}{tetlist(1)};            % retrieve ripples struct

times = r.timerange(1):0.001:r.timerange(end);              
nrip = zeros(size(times));

for t = 1:length(tetlist)                                   % iterate through trigger tetrodes
    tmprip = ripples{index(1,1)}{index(1,2)}{tetlist(t)};
    
    if (minthresh == 0)                                     % (if take any unextracted ripple)
        % get all the times
        rtimes = [tmprip.starttime tmprip.endtime];
    else
        % get the indeces for the ripples with energy above minthreshold
        rvalid = find(tmprip.maxthresh > minthresh);
        rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
    end
    
    % create another parallel vector (i.e. nrtimes) with bordering times for zeros
    nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];            % [ <starttimes>-.00001 <endttimes>+.00001 ] 
    rtimes = reshape(rtimes', length(rtimes(:)), 1);                        % [ <start1> ; <end1> ; <start2> ; <end2> ; .... ]
    rtimes(:,2) = 1;                                                        % adds a column of 1s
    nrtimes = [r.timerange(1) ; reshape(nrtimes', ...                       % nrtimes has bordering times
        length(nrtimes(:)), 1) ; r.timerange(2)];
    nrtimes(:,2) = 0;
    
    % create a new list with all of the times in it
    tlist = sortrows([rtimes ; nrtimes]);                                   % paired 1s and 0s
    
    % use interp to create a set of ones and zeros for each time
    % and add to nrip to get a cumulative count of the number of
    % ripples per timestep
    try
        nrip = nrip + interp1(tlist(:,1), tlist(:,2), times, 'nearest');      % nearest-neighbor interpolation of ripple times
    catch
        keyboard
    end
end

rip.times = times;  %same sampling as pos, spikes, etc
clear times;
rip.nripples = nrip; %number of ripples on each tetrode

%apply excludetimes
includetimes = ~isExcluded(rip.times, excludeperiods); %list of ones and zeros sampled every millisecond, ones = included, zeros = excluded
if size(rip.nripples) ~= size(includetimes)
    includetimes = includetimes';
end
includerips = rip.nripples .* includetimes;

%make a list of ripple start times
a = zeros(1, length(includerips));
a(find(includerips >= 1)) = 1;
ripstart = rip.times(diff(a)==1);

%Calculate the xcorrelation between spikes and ripple starttimes.
cells = unique(index(:,5:6),'rows');
for i = 1:length(cells(:,1))
    cell = cells(i,:);
    spiketimes = spikes{index(1,1)}{index(1,2)}{cell(1,1)}{cell(1,2)}.data(:,1);
    if ~isempty(ripstart)
    xcorrstruct = spikexcorr(spiketimes, ripstart, binsize, window(2));
    c1vsc2(i,:) = xcorrstruct.c1vsc2/sqrt(xcorrstruct.nspikes1*xcorrstruct.nspikes2);   % final output is normalize cross-correlogram
    end
end

if exist('c1vsc2')
out.c1vsc2 = c1vsc2;
out.time = xcorrstruct.time;
end
end