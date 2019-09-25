function out = DFAsj_HPexpt_getriprate(index, excludetimes, ripples, tetinfo, pos, varargin)
% out = DFAsj_getripalignspiking(spike_index, excludeperiods, spikes, ripples, tetinfo, options)

% Called from DFSsj_HPexpt_getriprate. Similar to DFAsj_getripalignspiking
% Use tetinfo and tetfilter passed in, or redefine here to get riptets
% Then use ripples to getriptimes. Can use inter-ripple-interval of 1 sec, and use a low-speed criterion.
% Can get epoch time period from ripples, or from pos


tetfilter = '';
excludetimes = [];
maxcell = 0;
minstd = 3;
lowsp_thrs = 5; %cm/sec
highsp_thrs = lowsp_thrs;
dospeed = 0;
doiri=0; % interripple-interval

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'tetfilter'
            tetfilter = varargin{option+1};
        case 'excludetimes'
            excludetimes = varargin{option+1};
        case 'minstd'
            minstd = varargin{option+1};
        case 'maxcell'
            maxcell = varargin{option+1};
        case 'dospeed'
            dospeed = varargin{option+1};
        case 'doiri'
            doiri = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

index;
day = index(1,1);
epoch = index(1,2);

% Get riptimes
% -------------
if isempty(tetfilter)
    riptimes = sj_getripples_tetinfo([day epoch], ripples, tetinfo, 'tetfilter', '(isequal($descrip, ''riptet''))','minstd',3);
else
    riptimes = sj_getripples_tetinfo([day epoch], ripples, tetinfo, 'tetfilter', tetfilter, 'minstd', 3);
end
% Can opt to have a cellcountthresh for each event as in getpopulationevents2 or  sj_HPexpt_ripalign_singlecell_getrip4
% Not using as of now

% Get triggers as rip starttimes separated by at least 1 sec
% ----------------------------------------------------------
rip_starttime = 1000*riptimes(:,1);  % in ms


if doiri
    % Find ripples separated by atleast a second
    % --------------------------------------------
    iri = diff(rip_starttime);
    keepidx = [1;find(iri>=1000)+1];
    rip_starttime = rip_starttime(keepidx);
end

% Implement speed criterion - Skip for now, or keep. Try both
% ----------------------------------------
if dospeed
    absvel = abs(pos{day}{epoch}.data(:,5)); % Can also use field 9
    postime = pos{day}{epoch}.data(:,1); % in secs
    pidx = lookup(rip_starttime,postime*1000);
    speed_atrip = absvel(pidx);
    lowsp_idx = find(speed_atrip <= lowsp_thrs);
    highsp_idx = find(speed_atrip > highsp_thrs);
    
    rip_starttime = rip_starttime(lowsp_idx);
end

% Get the time period of the epoch time
% --------------------------------------
postime = pos{day}{epoch}.data(:,1); % in secs
totaltime = postime(end)-postime(1); % in secs
riprate = length(rip_starttime)./totaltime; % Hz, ripples/ sec


% Output
% ------
out.index = index; % day epoch
out.riprate = riprate;
out.rip_starttime = rip_starttime;
out.totaltime = totaltime;

