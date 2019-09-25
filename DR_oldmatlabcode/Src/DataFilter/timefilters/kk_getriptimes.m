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
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'cellfilter'
            cellfilter = varargin{option+1};
        case 'minenergy'
            minenergy = varargin{option+1};
        case 'minthresh'
            minenergy = varargin{option+1};
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
        out{epochs(i,1)}{epochs(i,2)}.time = times;
        clear times;
        out{epochs(i,1)}{epochs(i,2)}.nripples = nrip;
        
    else    % empty tetlist this epoch
        disp(sprintf('getriptimes -- no valid tetrodes: day %d epoch %d',epochs(i,1),epochs(i,2)))
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
