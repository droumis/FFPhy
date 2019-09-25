function [out] = getgammatimes(animaldir,animalprefix, epochs, tetlist, type, varargin)
% out = getgammatimes(animaldir,animalprefix,epochs, tetlist, options)
%
%     animaldir and animal prefix are strings indicating the base director for
%     the animal's data and the prefix for the data files
%
%     epochs is an Nx2 list of days and epochs
%
%     tetlist is a list of tetrodes to use or an empty matrix if the
%     'cellfilter' option is used.
%     
%     type identifies whether low or high gamma times should be used.
%
% options are
%	'cellfilter', 'cellfilterstring'
%		     specifies a cell filter to select the tetrodes to use for
%		     ripple filtering
%   'tetlist', 'tetfilterstring'
%            specifies a tetfilter to select the tetrodes to use for gamma
%            filtering
%	'minthresh', minthresh 
%		     specifies a minimum threshold in stdev units for a valid 
%			gamma event  (default 0)
%   'inclusive', 0 or 1
%           if set to 1, will link together gammatimes that are nearby.
%           (default 0)
% Produces a cell structure with a time field and an ngammas field which
% indicates the number of electrodes with a gamma event at each time point
%
% Examples:
% getgammatimes('/data/name/Fre', 'fre', epochs, 1,'high')
% getgammatimes('/data/name/Fre', 'fre', epochs, [], 'cellfilter', '(isequal($area, ''CA1''))')

% assign the options
cellfilter = '';
tetfilter = '';
minthresh = 0;
inclusive = 0;
min_separation = 1;

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'cellfilter'
            cellfilter = varargin{option+1};
        case 'tetfilter'
            tetfilter = varargin{option+1};
        case 'minthresh'
            minthresh = varargin{option+1};
        case 'inclusive'
            inclusive = varargin{option+1};
        case 'min_separation'
            min_separation = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

%check to see if a cell filter is specified
if ~isempty(cellfilter)
    % this will cause us to ignore tetlist
    cellinfo = loaddatastruct(animaldir, animalprefix, 'cellinfo');
end

%check to see if a tetrode filter is specified
if ~isempty(tetfilter)
    tetinfo = loaddatastruct(animaldir,animalprefix,'tetinfo');
end

loaddays = unique(epochs(:,1));
if isequal('high',type)
    gam = loaddatastruct(animaldir, animalprefix, 'gammah', loaddays);
elseif isequal('low',type)
    gam = loaddatastruct(animaldir, animalprefix, 'gammal', loaddays);
end

for i = 1:size(epochs,1)
    % if cellfilter is set, apply it to this day and epoch
    if ~isempty(cellfilter)
        tetlist =  evaluatefilter(cellinfo{epochs(i,1)}{epochs(i,2)}, ...
            cellfilter); 
        % get rid of the cell indeces and extract only the tetrode numbers 
        tetlist = unique(tetlist(:,1))';
    end
    % if tetfilter is set, apply it to this day and epoch
    if ~isempty(tetfilter)
        tetlist = evaluatefilter(tetinfo{epochs(i,1)}{epochs(i,2)},tetfilter);
    end
    if (~isempty(tetlist))
	% go through the tetlist and construct an an array where each element 
	% represents the number of active tetrodes for each 1 ms timestep.
	try
	    g = gam{epochs(i,1)}{epochs(i,2)}{tetlist(1)};
	catch
	    keyboard
	end
	times = g.timerange(1):0.001:g.timerange(end);
	ngam = zeros(size(times));
	for t = 1:length(tetlist)
	    tmpgam = gam{epochs(i,1)}{epochs(i,2)}{tetlist(t)};
	    % apply the minthresh threhsold
        if ~isempty(tmpgam.startind)
    	    gvalid = find(tmpgam.maxthresh > minthresh);
            if ~isempty(gvalid)
                if inclusive
                    tmp_start = tmpgam.starttime(gvalid);
                    tmp_end = tmpgam.endtime(gvalid);
                    tmp_index = nan(length(tmp_start),2);
                    tmp_index(1,1) = tmp_start(1);
                    count = 1;
                    for a = 2:length(tmp_start)
                        if (tmp_start(a) - tmp_index(count,1)) > min_separation
                           tmp_index(count,2) = tmp_end(a-1);
                           count = count+1;
                           tmp_index(count,1) = tmp_start(a);
                        end
                        if a == length(tmp_start)
                            tmp_index(count,2) = tmp_end(end);
                        end
                    end
                    tmp_start = tmp_index(~isnan(tmp_index(:,1)),1);
                    tmp_end = tmp_index(~isnan(tmp_index(:,2)),2);
                    gtimes = [tmp_start tmp_end];
                    clear tmp_index tmp_start tmp_end
                else
                    gtimes = [tmpgam.starttime(gvalid) tmpgam.endtime(gvalid)];
                end
                % create another parallel vector with bordering times for zeros
                ngtimes = [(gtimes(:,1) - 0.00001) (gtimes(:,2) + 0.00001)];
                gtimes = reshape(gtimes', length(gtimes(:)), 1); 
                gtimes(:,2) = 1;
                ngtimes = [g.timerange(1) ; reshape(ngtimes', ...
                    length(ngtimes(:)), 1) ; g.timerange(2)];
                ngtimes(:,2) = 0;
                % create a new list with all of the times in it
                tlist = sortrows([gtimes ; ngtimes]);
                %Get rid of any repeated times
                [b,m,n] = unique(tlist(:,1));
                tlist = [tlist(m,1) tlist(m,2)];
                clear b m n
                try
                    ngam = ngam + interp1(tlist(:,1),tlist(:,2),times,'nearest');
                catch
                keyboard
                end
            end
        end
	end
	out{epochs(i,1)}{epochs(i,2)}.time = times;
	clear times;
	out{epochs(i,1)}{epochs(i,2)}.ngamma = ngam;
    end
end
