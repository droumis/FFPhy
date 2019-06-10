function [out] = dfa_riptrigspiking(index, excludeperiods, spikes, ...
    eventscons, pos, task, varargin)

% MS adapted 2016 from dfakk_geteventtriggeredspiking
% This function finds spiking triggered to LFP events such as ripples.

disp(sprintf('%d %d %d %d',index))

% Note that event times are inherited from kk_getconstimes
    % (thus kk_getconstimes is where you should specify minthresh (i.e. 3
    % SD or 7 SD etc) - parameters can be entered in the DFScript

%
%   index [day epoch tetrode cell tetrode cell]
%
%%% Options are
%	'minthresh'-- specifies the minimum threshold of a valid event
%   'window'   -- 1x2 vector specifies the window before and after each included event.
%                   Default is 500 mseconds before and after
%                   event start time.
%   'TF'       -- tetfilter, use to filter for tetrodes for which to find
%                    consensus events
%   'binsize'  -- histc binsize for spike rasters, 1 ms
%   'frbinsize'-- binsize to calculate the firing rate, default 10 ms
%   'welldist' -- to specify that events should only be included if they
%                 occur a certain distance from the wells
%   The remaining options are for getconstimes.

%   out = out.out       An R x C sized matrix where each entry is the number of
%                       spikes cell C fired during event R
%         out.index     [D E T C], gives the identity of the cells for each
%                       column in out.out
%         out.time      vector of binedges used for histc
%         out.psth      Array of 1s and 0s, histc count of spikes around
%                       ripples for each ripple in 1 ms bins
%         out.frpsth    Array of histc spikes counts in larger bins around
%                       each ripple

% default options
minvelocity = 0;  
welldist = [];
window = [0.5 0.5]; % in sec
binsize = 0.001; % 1 ms for rasters
frbinsize= 0.01; % 10 ms for population FR plotting
TF = 1;
% % Set options
if (~isempty(varargin))
    assign(varargin{:});
end
% for option = 1:2:length(varargin)-1
%     
%     switch varargin{option}
%         case 'window'
%             window = varargin{option+1};
%         case 'TF'                % identify tetrodes for detecting events if not using consensus
%             TF = varargin{option+1};   
%         case 'binsize'
%             binsize = varargin{option+1};
%         case 'frbinsize'
%             frbinsize = varargin{option+1};
%         case 'minthresh'
%             minthresh = varargin{option+1};            
%         case 'consensus_numtets'
%             consensus_numtets = varargin{option+1};
%         case 'maxvelocity'
%             maxvelocity = varargin{option+1};
%         case 'minvelocity',
%             minvelocity = varargin{option+1};
%         case 'welldist',
%             welldist = varargin{option+1};            
%         otherwise
%             error(['Option ''', varargin{option}, ''' not defined']);
%     end
% end

emptyoutput_flag = 0;

day = index(1);
epoch = index(2);

% animaldir,animalprefix

% First receive valid event periods for the epoch's day.

ec = eventscons{day}{epoch}{TF};

% if the eventcons is empty or not found..
if ~exist('ec','var')
    ec = [];
    includetimes = [];   
else
    try 
        times = ec.eegtimesvec_ref; %try to use the lfp times from the ec struct if it exists
    catch
        times = ec.timerange(1):1/ec.samprate:ec.timerange(end);
    end
    %obtain includetimes
    includetimes = ~isExcluded(times, excludeperiods);     % list of ones and zeros sampled every millisecond, ones = included, zeros = excluded
        includetimes = includetimes(:);
end

% empty output checking, 2 scenarios: %%%%%%%%%%%%%%%%%%%%%%%%
try
if sum(includetimes)==0 || isempty(ec) || isempty(ec.eventname)
    disp(sprintf('dfa: no includetimes d%de%d received from kk_getconstimes',day,epoch));
  % also, on some odd old animals (dudley..) cellinfo has some entries while
    % the spikes var does not have the corresponding entry.. this is to
    % continue past these
    emptyoutput_flag = 1;
elseif index(4) > length(spikes{index(1)}{index(2)}{index(3)}) || ...
       isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
    disp('cellinfo has an entry that spikes does not -- ignoring..')
    emptyoutput_flag = 1;
end
catch
    fprintf('no events. skipping\n');
    emptyoutput_flag = 1;
end
if emptyoutput_flag
    out.index = index;
%     out.type = '';
    out.epoch_type = '';
    out.epoch_environment = '';
    out.time = [];
    out.frtime = [];
    out.psth = [];
    out.frhist = [];
    out.instantFR = [];
    out.nospikes = [];
    out.noevents = [];
    out.eventtags = [];
    out.posteventmatrix = [];
    out.eventduration = [];
%     out.sleepc_nospikes = nan;
%     out.sleepc_totalduration = nan;
%     out.sleepc_time_immobile = nan;
%     out.sleepc_velocity_thresh = nan;      
    return   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(ec.tetlist) < consensus_numtets
    error('something is wrong-- your eventcons data should have at least minimum number of consensus tets')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate psths ---------

% first reconstitute the cons event (originally from kk_getconstimes)
eventstart = times([0 diff(includetimes')]==1);
eventend = times(diff(includetimes')==-1);

% throw out last event if interrupted by end of epoch
if (eventend(end)-eventstart(end) < 0)
    eventstart = eventstart(1:(end-1));
end
% throw out first event if occurred before beginning of epoch
if (eventend(1)-eventstart(1) < 0)
   eventend = eventend(2:end); 
end
% throw out any event that begins less than window(2) (i.e. 0.5 seconds) from end of epoch
while (eventend(end) > times(end)-window(2)) 
    eventend(end) = [];
    eventstart(end) = [];
end

eventstartend = [eventstart' eventend'];

% choose event time to be the START of the event period
eventtimes = eventstartend(:,1);

% (optional) if specified, filter for ripples that occur > welldist from well
% this is in reference to the start time of the ripple
if ~isempty(welldist)
     o = kk_getwelltimes(animaldir,animalprefix,[day epoch],welldist);  % returns nan if not a linear run epoch w/ linpos data
     if ~isnan(o{day}{epoch}.time)
         nonwellperiods = vec2list(~o{day}{epoch}.nearwell,o{day}{epoch}.time);
         eventtimes = eventtimes(logical(isExcluded(eventtimes,nonwellperiods)));
         eventstartend = eventstartend(logical(isExcluded(eventstartend(:,1),nonwellperiods)), :);
     end
end

% Now obtain posteventmatrix (0 and 1s) of consensus times
    % basically rehashing kk_getconstimes, except taking all cons times
    % (instead of excluding chained / exclusion2-overlapping ones) since we want to see them later
postevent_numbins = window(2)/binsize;
posteventmatrix = nan(length(eventtimes),postevent_numbins);
    % First obtain consensus events
        % before doing so, we need to filter by size + max velocity + min velocity
    eclist = [ec.starttime ec.endtime];
        % size
    eclist = eclist((ec.maxthresh > minthresh),:);
        % in specified velocity range
    posentry = pos{day}{epoch};
    posinds = lookup(eclist(:,1),posentry.data(:,1));  % velocity at starttime
    if size(posentry.data,2) > 5
        goodinds = (posentry.data(posinds,9) < maxvelocity) & (posentry.data(posinds,9) >= minvelocity);
    else
        goodinds = (posentry.data(posinds,5) < maxvelocity)  & (posentry.data(posinds,5) >= minvelocity);
        disp(sprintf('.. using column 5 velocity..'))
    end
    eclist = eclist(goodinds,:); 
    ecvec = list2vec(eclist,times)';
    % Now transcribe into posteventmatrix
    for rr = 1:length(eventtimes)
        ind = lookup(eventtimes(rr),times);
        posteventmatrix(rr,:) = ecvec(ind:(ind+postevent_numbins-1));
    end

   
%Calculate PSTHs.

% histc bins, where center bin is centered time 0
time = (-window(1)-0.5*binsize):binsize:(window(2)+0.5*binsize);
frtime = (-window(1)-0.5*frbinsize):frbinsize:(window(2)+0.5*frbinsize);
% retrieve spikes
if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
    spiketimes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1);
else
    spiketimes = [];
end

% now calculate histograms (while also creating epoch and time "tag" for events)
if ~isempty(eventtimes)
    
    psth = nan(length(eventtimes),length(time));
    frhist = nan(length(eventtimes),length(frtime));
    eventtags = nan(length(eventtimes),2);          %  [epoch  ,  eventtime]       % useful later on for plotting and more
    instantFR = nan(length(eventtimes),length(time));
    
    % iterate through each event and calculate psth
    for r=1:length(eventtimes)
        % psth and instantaneous firing rate 
        onehist = histc(spiketimes , eventtimes(r) + time);
        twohist = histc(spiketimes , eventtimes(r) + frtime);
        ifr = instantfr(spiketimes, [eventtimes(r) + time]);
        % get concatenation in right dimension
        if isempty(onehist)
            psth(r,:) = zeros(size(time));
            eventtags(r,:) = [epoch eventtimes(r)];
            instantFR(r,:) = zeros(size(time));
            frhist(r,:) = zeros(size(frtime));

        elseif size(onehist,2) == length(time)
            psth(r,:) = onehist;
            eventtags(r,:) = [epoch eventtimes(r)];
            instantFR(r,:) = ifr;
            frhist(r,:) = twohist;
        else
            psth(r,:) = onehist';
            eventtags(r,:) = [epoch eventtimes(r)];
            instantFR(r,:) = ifr;
            frhist(r,:) = twohist';
        end
    end
    
    % if available, add sleep classification to eventtags  "sleepc" event
    %    (0: not included in a sleep period, while 1: is included)
%     if (length(sleep{day}) >= epoch)  &&  ~isempty(sleep{day}{epoch})
%         sleepincludetimes  =   [  sleep{day}{epoch}.starttime     sleep{day}{epoch}.endtime  ];
%         eventtags = [  eventtags   isExcluded(eventtimes,sleepincludetimes)  ];
%     else
%         eventtags = [  eventtags   zeros(size(eventtags,1),1)   ];
%     end
end

%%% outputs %%%%%%%%%
out.index = index;
% out.type = cellinfo{index(1)}{index(2)}{index(3)}{index(4)}.type;
out.epoch_type = task{index(1)}{index(2)}.type;
if isfield(task{index(1)}{index(2)},'environment')
    out.epoch_environment = task{index(1)}{index(2)}.environment;
else
    out.epoch_environment = '';
end
out.time = time;
out.eventstartend = ;
out.frtime = frtime;
% out.psth = sparse(psth);
out.psth = psth;
out.frhist = frhist;
out.instantFR = instantFR;
out.nospikes = length(spiketimes);      % number of spikes in whole epoch
out.noevents = length(eventtimes);      % this is the number of events this epoch (kk_get<event>times should
                                       % exclude events that come too soon after first event via 'exclusion'
                                       % option)
out.eventtags = eventtags;             % [ <epoch no>  <time of event occurrence>   <sleepc event, == 1 if it is> ]
out.posteventmatrix = posteventmatrix;             % [ <epoch no>  <time of event occurrence>   <sleepc event, == 1 if it is> ]
out.eventduration =  eventend'-eventstart';
% if strcmp(out.epoch_type,'sleep') 
%     out.sleepc_nospikes = sum(isExcluded(spiketimes,sleepincludetimes));     %   # of spikes that fall into "sleep-classification" immobility periods
%     if ~isempty(sleepincludetimes)
%         out.sleepc_totalduration = sum(sleepincludetimes(:,2)-sleepincludetimes(:,1));
%     else
%         out.sleepc_totalduration = 0;
%     end
%     out.sleepc_time_immobile = sleep{day}{epoch}.time_immobile;
%     out.sleepc_velocity_thresh = sleep{day}{epoch}.velocity_thresh;
% else
%     out.sleepc_nospikes = nan;
%     out.sleepc_totalduration = nan;
%     out.sleepc_time_immobile = nan;
%     out.sleepc_velocity_thresh = nan;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end