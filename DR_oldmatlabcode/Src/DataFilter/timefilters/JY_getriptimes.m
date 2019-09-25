function [out] = JY_getriptimes(animaldir,animalprefix, epochs, tetlist, minrip,varargin)
% out = getriptimes(animaldir,animalprefix,epochs, tetlist, options)
%  
%     animaldir and animal prefix are strings indicating the base director for
%     the animal's data and the prefix for the data files
%
%     epochs is an Nx2 list of days and epochs
%
%     tetlist is a list of tetrodes to use or an empty matrix if the
%     'cellfilter' option is used.
% minrip -  more than this number of tetrodes need to have ripples for the
%           event to be counted
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
minthresh = 3; % 3SD
% minimum number of simultaneous ripples
minrip=minrip;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'cellfilter'
            cellfilter = varargin{option+1};
        case 'minenergy'
            minenergy = varargin{option+1};
        case 'minthresh'
            minthresh = varargin{option+1};
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
        % get rid of the cell indeces and extract only the tetrode numbers
        tetlist = unique(tetlist(:,1))';
    end
    if (~isempty(tetlist))
        % go through the tetlist and construct an an array where each element
        % represents the number of active tetrodes for each 1 ms timestep.
        try
            r = rip{epochs(i,1)}{epochs(i,2)}{tetlist(1)};
        catch
            display('error');
            keyboard
        end
        times = r.timerange(1):0.001:r.timerange(end);
        nrip = zeros(size(times));
        for t = 1:length(tetlist)
            tmprip = rip{epochs(i,1)}{epochs(i,2)}{tetlist(t)};
            % apply the minthresh threhsold
            rvalid = find(tmprip.maxthresh > minthresh);
            rtimes = [tmprip.starttime(rvalid) tmprip.endtime(rvalid)];
            % create another parallel vector with bordering times for zeros
            nrtimes = [(rtimes(:,1) - 0.00001) (rtimes(:,2) + 0.00001)];
            rtimes = reshape(rtimes', length(rtimes(:)), 1);
            rtimes(:,2) = 1;
            nrtimes = [r.timerange(1) ; reshape(nrtimes', ...
                length(nrtimes(:)), 1) ; r.timerange(2)+0.00002];
            nrtimes(:,2) = 0;
            % create a new list with all of the times in it
            tlist = sortrows([rtimes ; nrtimes]);
            % use interp to create a set of ones and zeros for each time
            % and add to nrip to get a cumulative count of the number of
            % ripples per timestep
            try
                nrip = nrip + interp1(tlist(:,1), tlist(:,2), times, 'nearest');
            catch
                % we get here if there are multiple tlist entries with the same
                % time, so we need to get rid of them
                
                keyboard
            end
            
            %find the start and end borders of each ripple
            %     inripple = (nrip >= minrip);
            %     startrippleind = find(diff(inripple) == 1)+1;
            %     endrippleind = find(diff(inripple) == -1)+1;
            %     ripplestdout = [];
            %     if ((length(startrippleind) > 1) & (length(endrippleind) > 1))
            %         if (endrippleind(1) < startrippleind(1))
            %             endrippleind = endrippleind(2:end);
            %         end
            %         if (endrippleind(end) < startrippleind(end))
            %             startrippleind = startrippleind(1:end-1);
            %         end
            %         startripple = times(startrippleind);
            %         endripple = times(endrippleind);
            % %         for i = 1:length(startripple)
            % %             ripplestdout(i,1) = max(ripplestd(startrippleind(i):endrippleind(i)));
            % %         end
            %
            %         out = [startripple(:) endripple(:)];
            %
            %         out = out( ((out(:,2)-out(:,1))>.050),:);
            %     else
            %         out = [];
            %         ripplestdout = [];
            %     end
            %
        end
        out{epochs(i,1)}{epochs(i,2)}.time = times;
        clear times;
        out{epochs(i,1)}{epochs(i,2)}.nripples = nrip;
        
        %end
        % 	out{epochs(i,1)}{epochs(i,2)}.time = times;
        % 	clear times;
        % 	out{epochs(i,1)}{epochs(i,2)}.nripples = nrip;
    end
end
