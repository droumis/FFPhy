function [out] = kk_getriptimes(animaldir,animalprefix, epochs, tetlist, varargin)
% out = getriptimes(animaldir,animalprefix,epochs, tetlist, options)

%%    modified by kk May 2013 to avoid end ripple time problems
%
%     animaldir and animal prefix are strings indicating the base director for
%     the animal's data and the prefix for the data files
%
%     epochs is an Nx2 list of days and epochs
%
%     tetlist is a list of tetrodes to use or an empty matrix if the
%     'cellfilter' option is used.
%
% options are
%	'cellfilter', 'cellfilterstring'
%		     specifies a cell filter to select the tetrodes to use for
%		     ripple filtering
%	'minthresh', minthresh 
%		     specifies a minimum threshold in stdev units for a valid 
%			ripple event  (default 0)
%
% Produces a cell structure with a time field and an nripples field which
% indicates the number of electrodes with a ripple at each time point
%
% Examples:
% getriptimes('/data/name/Fre', 'fre', epochs, 1)
% getriptimes('/data/name/Fre', 'fre', epochs, [], 'cellfilter', '(isequal($area, ''CA1''))')

% assign the options
cellfilter = '';
minenergy = 0;
minthresh = 0;
exclusion_dur = 0;     % duration after ripple within which start of any subsequent ripple is eliminated
exclusion_nrip = 0;    % number of tetrodes that must have ripple for exclusion to kick in

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'cellfilter'
            cellfilter = varargin{option+1};
        case 'minenergy'
            minenergy = varargin{option+1};
        case 'minthresh'
            minthresh = varargin{option+1};
        case 'exclusion'                 % kk implemented 8.1.13, excludes ripples 
            exclusion_dur = varargin{option+1};     % duration after ripple within which start of any subsequent ripple is eliminated
            exclusion_nrip = varargin{option+2};    % number of tetrodes that must have ripple for exclusion to kick in
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%check to see if a cell filter is specified
if (~isempty(cellfilter))
    % this will cause us to ignore tetlist
    cellinfo = loaddatastruct(animaldir, animalprefix, 'cellinfo');
end


loaddays = unique(epochs(:,1));
rip = loaddatastruct(animaldir, animalprefix, 'ripples', loaddays);
for i = 1:size(epochs,1)
    % if cellfilter is set, apply it to this day and epoch
    if (~isempty(cellfilter))
	tetlist =  evaluatefilter(cellinfo{epochs(i,1)}{epochs(i,2)}, ...
				cellfilter); 
    end
    if (~isempty(tetlist))
    % get rid of the cell indeces and extract only the tetrode numbers 
        tetlist = unique(tetlist(:,1))';
	% go through the tetlist and construct an an array where each element 
	% represents the number of active tetrodes for each 1 ms timestep.
	try
	    r = rip{epochs(i,1)}{epochs(i,2)}{tetlist(1)};
	catch
	    keyboard
    end
	times = r.timerange(1):0.001:r.timerange(end);
	nrip = zeros(size(times));
	for t = 1:length(tetlist)
	    tmprip = rip{epochs(i,1)}{epochs(i,2)}{tetlist(t)};
	    % apply the minthresh threhsold
	    rvalid = find(tmprip.maxthresh > minthresh);
        
	    rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
        
        % possibly no ripples.. if so skip (kk added 7.29.13)
        if isempty(rtimes)
            continue
        end
        
        %check for possible last ripple that extends past times vector
        if times(end)-rtimes(end,2) < 0
            rtimes(end,:) = [];
            disp(sprintf('excluded last ripple, day %d epoch %d',epochs(i,1),epochs(i,2)))
        end
        
	    % create another parallel vector (nrtimes) with bordering times for zeros
	    nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];     % borders each period with .01 ms
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
    
    % OUTPUT.
        % If applying exclusion window, then cull offending ripples  kk 8.1.13
        % Note that for clarity, we end up "flattening" surviving nrip to
        % exclusion_nrip.
    out{epochs(i,1)}{epochs(i,2)}.time = times;
    
    if exclusion_dur == 0
    
        clear times
        out{epochs(i,1)}{epochs(i,2)}.nripples = nrip;
    
    elseif exclusion_dur > 0
        
        nrip2 = nrip >= exclusion_nrip;
        
        if sum(nrip2) == 0
            disp(sprintf('kk_getriptimes: day %d epoch %d no periods satisfy exclusion_nrip',epochs(i,1),epochs(i,2)))
            out{epochs(i,1)}{epochs(i,2)}.nripples = zeros(size(times));
            clear times
            continue
        end
        
        %convert to [starttime endtime] 
        riplist = getExcludePeriods(times',~nrip2);
        % now iterate through each ripple period
        for rr=fliplr(2:size(riplist,1))
            timesinceprevrip = riplist(rr,1) - riplist(rr-1,2);
            if timesinceprevrip < exclusion_dur
                riplist(rr,:) = [ ];
                disp(sprintf('kk_getriptimes: chain-excluded ripple: d %d e %d (%d s)',epochs(i,1),epochs(i,2),timesinceprevrip))
            end
        end
        % convert back to 1 ms bin vector form (simple this time, since
        % times perfectly match)
        
        nrip = zeros(size(times));          % initialize
        
        % install 1s
        for k=1:size(riplist,1)
            startindex = find(riplist(k,1)==times);
            endindex = find(riplist(k,2)==times);
            nrip(startindex:endindex) = 1;
        end
        
        nrip = exclusion_nrip * nrip;    % assign threshold value
        
        out{epochs(i,1)}{epochs(i,2)}.nripples = nrip;
        clear times;
        
    end
    
    
        
    
    
    else    % empty tetlist this epoch
        disp(sprintf('kk_getriptimes: no valid tetrodes: day %d epoch %d',epochs(i,1),epochs(i,2)))
        % retrieve some empty times structure for that epoch
        for tet=1:7
            if ~isempty(rip{epochs(i,1)}{epochs(i,2)}{tet})
                r = rip{epochs(i,1)}{epochs(i,2)}{tet};
                times = r.timerange(1):0.001:r.timerange(end);
                continue
            end
        end
        out{epochs(i,1)}{epochs(i,2)}.time = times;
        out{epochs(i,1)}{epochs(i,2)}.nripples = zeros(size(times));
        clear times;
    end
    
    
    
end
