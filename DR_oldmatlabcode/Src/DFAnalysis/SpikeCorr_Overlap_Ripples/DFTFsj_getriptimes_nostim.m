function [out] = DFTFsj_getriptimes_nostim(animaldir,animalprefix, epochs, varargin)
% out = getriptimes(animaldir,animalprefix,epochs, tetlist, options)
% Shantanu: Like getriptimes, but uses DIO to set timewin around stimtimes
% to 0. Needs DIO
%
%     animaldir and animal prefix are strings indicating the base director for
%     the animal's data and the prefix for the data files
%
%     epochs is an Nx2 list of days and epochs
%
%     tetrodes specified as tetrode list or tetrode filter or cellfilter
%     'cellfilter'or 'tetfilter' option is used.
%
% options are
%	'cellfilter', 'cellfilterstring'
%		     specifies a cell filter to select the tetrodes to use for
%		     ripple filtering
%   'tetfilter', 'tetfilterstring'
%            specifies a tetfilter to select the tetrodes to use for ripple
%            filtering
%	'minthresh', minthresh 
%		     specifies a minimum threshold in stdev units for a valid 
%			ripple event  (default 0)
% Produces a cell structure with a time field and an nripples field which
% indicates the number of electrodes with a ripple at each time point
%
% Examples:
% getriptimes('/data/name/Fre', 'fre', epochs, 1)
% getriptimes('/data/name/Fre', 'fre', epochs, [], 'cellfilter', '(isequal($area, ''CA1''))')

% assign the options
cellfilter = '';
tetfilter = '';
minthresh = 0;
tetrodes = [];
timewin = 0.1; % Time Window around stimtime to exclude  100ms - 0.1sec 
timewin1 = []; % Time Window before stimtime to exclude  100ms - 0.1sec
timewin2 = []; % Time Window after stimtime to exclude  100ms - 0.1sec
excludeperiods = [];
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'tetrodes'
            tetrodes = varargin{option+1};
        case 'cellfilter'
            cellfilter = varargin{option+1};
        case 'tetfilter'
            tetfilter = varargin{option+1};
        case 'minthresh'
            minthresh = varargin{option+1};
        case 'timewin'
            timewin = varargin{option+1};
        case 'timewin1'
            timewin1 = varargin{option+1};   
        case 'timewin2'
            timewin2 = varargin{option+1};
        case 'excludeperiods'
            excludeperiods = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

if isempty(timewin1)
    timewin1 = timewin;
end
if isempty(timewin2)
    timewin2 = timewin;
end


% Tetrode selection
if ~isempty(tetrodes)
    tetlist = tetrodes;
end
%check to see if a cell filter is specified
if ~isempty(cellfilter)
    % this will cause us to ignore tetlist if specified
    cellinfo = loaddatastruct(animaldir, animalprefix, 'cellinfo');
end
%check to see if a tetrode filter is specified
if ~isempty(tetfilter)
    % this will cause us to ignore tetlist if specified
    tetinfo = loaddatastruct(animaldir,animalprefix,'tetinfo');
end

loaddays = unique(epochs(:,1));
rip = loaddatastruct(animaldir, animalprefix, 'ripples', loaddays);
DIO = loaddatastruct(animaldir, animalprefix, 'DIO', loaddays);

for i = 1:size(epochs,1)
    
    % Get stimtimes
    %--------------
    try
        stim = DIO{epochs(i,1)}{epochs(i,2)}{16};
        if isempty(stim)
            stim = DIO{day}{epoch}{15};
        end
        stim_starttime = stim.pulsetimes(:,1)./10000; %sec (divide by 10 is ms)
        stim_endtime = stim.pulsetimes(:,2)./10000; %sec
    catch
        % Should set stim to empty. But keep keyboard for now for debug
        keyboard
        stim_starttime = [];
    end
    
    
    % Get ripples for timerange
    %------------
    % if cellfilter is set, apply it to this day and epoch - Ignore tetlist
    if ~isempty(cellfilter)
        tetlist =  evaluatefilter(cellinfo{epochs(i,1)}{epochs(i,2)}, ...
            cellfilter);
        % get rid of the cell indeces and extract only the tetrode numbers
        tetlist = unique(tetlist(:,1))';
    end
    % if tetfilter is set, apply it to this day and epoch - Ignore tetlist
    if ~isempty(tetfilter)
        tetlist = evaluatefilter(tetinfo{epochs(i,1)}{epochs(i,2)},tetfilter);
    end
    if (~isempty(tetlist))
        % go through the tetlist and construct an an array where each element
        % represents the number of active tetrodes for each 1 ms timestep.
        try
            r = rip{epochs(i,1)}{epochs(i,2)}{tetlist(1)}; % For timerange
        catch
            keyboard
        end
        
        times = r.timerange(1):0.001:r.timerange(end);
        nrip = zeros(size(times));
        
        % Set nrip based on ripples on tetrodes in telist
        for t = 1:length(tetlist)
            tmprip = rip{epochs(i,1)}{epochs(i,2)}{tetlist(t)};
            % apply the minthresh threhsold
            if ~isempty(tmprip.startind)
                rvalid = find(tmprip.maxthresh > minthresh);
                if ~isempty(rvalid)
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
            end
        end
        
        % Set timewin around stimtime to 0 in nrip field 
        %-------------------------------------------
        % 1st method: Direct
        if (~isempty(stim_starttime))
            rtimes = [stim_starttime-timewin1  stim_starttime+timewin2];
            stidxs = lookup(rtimes(:,1), times);
            endidxs = lookup(rtimes(:,2), times);
            for s=1:length(stidxs),
                nrip(stidxs(s):endidxs(s))=0;
            end
        end
        
        % Update out
        out{epochs(i,1)}{epochs(i,2)}.time = times;
        clear times;
        out{epochs(i,1)}{epochs(i,2)}.nripples = nrip;
    end
end







