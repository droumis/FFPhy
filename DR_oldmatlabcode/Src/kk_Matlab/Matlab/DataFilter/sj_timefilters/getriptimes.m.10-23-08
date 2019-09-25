function [out] = getriptimes(animaldir,animalprefix, epochs, tetlist, varargin)
% out = getriptimes(animaldir,animalprefix,epochs, tetlist, options)
%
%     animaldir and animal prefix are strings indicating the base directory for
%     the animal's data and the prefix for the data files
%
%     epochs is an Nx2 list of days and epochs
%
%     tetlist is a list of tetrodes to use or an empty matrix if the
%     'cellfilter option is used.
%
% options are
%	'cellfilter', 'cellfilterstring'
%		     specifies a cell filter to select the tetrodes to use for
%		     ripple filtering
%	'minenergy', E
%		     specifies the minimum energy of a valid ripple event
%   'nsimul', minimum number of tetrodes that simultaneously recorded a ripple to
%           include it in analysis
%   'maxcell', = 1: find ripples on only one tetrode of those include in tetlist or
%           cellfilter, the tetrode with maximum number of cells
%
%
% Produces a cell structure with a time field and an nripples field which
% indicates the number of electrodes with a ripple at each time point
%
% Examples:
% getriptimes('/data/name/Fre', 'fre', epochs, 1)
%	sets the exclude list for the specified epochs to exclude non-ripple
%	events detected on tetrode 1
% getriptimes('/data/name/Fre', 'fre', epochs, [], 'nsimul', 2, 'cellfilter', '(isequal($area, ''CA1''))')
%	sets the exclude list for the specified epochs so that only times where
%	two or more CA1 tetrodes have a ripple are included.

% assign the options
cellfilter = [];
minenergy = 0;
maxcell = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'cellfilter'
            cellfilter = varargin{option+1};
        case 'minenergy'
            minenergy = varargin{option+1};
        case 'maxcell'
            maxcell = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end
keyboard

%check to see if a cell filter is specified
if (~isempty(cellfilter))
    % this will cause us to ignore tetlist
    cellinfo = loaddatastruct(animaldir, animalprefix, 'cellinfo');
end

%find tetrode with max number of cells per day
if maxcell == 1
    tetlistdays = [];
    cellinfo = loaddatastruct(animaldir, animalprefix, 'cellinfo');
    task = loaddatastruct(animaldir, animalprefix, 'task');
    a = sortrows(epochs,1);
    days = unique(a(:,1));
    for i = 1:length(days)  %find the tetrode with most cells on each day
        d = days(i);
        numcells = [0 0];
        tet =  gettetmaxcell(cellinfo, task, d, cellfilter, 2); %for all run epochs
        tetlistdays = [tetlistdays; d tet];
    end
end


loaddays = unique(epochs(:,1));
rip = loaddatastruct(animaldir, animalprefix, 'ripples', loaddays);
for i = 1:size(epochs,1)
    % if cellfilter is set, apply it to this day and epoch
    if (~isempty(cellfilter) & maxcell==0)
        tetlist =  evaluatefilter(cellinfo{epochs(i,1)}{epochs(i,2)}, ...
            cellfilter);
        % get rid of the cell indeces and extract only the tetrode numbers
        tetlist = unique(tetlist(:,1))';
    elseif maxcell == 1
        tetlist = tetlistdays(tetlistdays(:,1)==epochs(i,1),2);
    end
    if (~isempty(tetlist))
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
        out{epochs(i,1)}{epochs(i,2)}.time = times';
        clear times;
        out{epochs(i,1)}{epochs(i,2)}.nripples = nrip';
    end
end
end