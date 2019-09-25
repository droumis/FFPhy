function [out] = kk_getgammaltimes(animaldir,animalprefix, epochs, tetlist, varargin)
% out = getriptimes(animaldir,animalprefix,epochs, tetlist, options)

%%    stupidly adapted from kk_getriptimes on 7.31.13
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
tetfilter = '';
minenergy = 0;
minthresh = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'cellfilter'
            cellfilter = varargin{option+1};
        case 'tetfilter'
            tetfilter = varargin{option+1};
        case 'minenergy'
            minenergy = varargin{option+1};
        case 'minthresh'
            minthresh = varargin{option+1};
        case 'ngammal_report'          % really used below just for the reporting of % epoch that qualifies
            ngammal_report = varargin{option+1};            
        case 'onephasetet'
            onephasetet = varargin{option+1};   % flag to restrict one phasetet a day (first one with most # of cells)            
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

% load cellinfo or tetinfo

if (~isempty(cellfilter))
    % this will cause us to ignore tetlist
    cellinfo = loaddatastruct(animaldir, animalprefix, 'cellinfo');
elseif (~isempty(tetfilter))
    % this will cause us to ignore tetlist
    tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo');
end

loaddays = unique(epochs(:,1));
rip = loaddatastruct(animaldir, animalprefix, 'gammal', loaddays);

    % if onephasetet specified, then obtain list of principal tetrodes to
    % use for each day (tetrode w/ most # cells)
    if onephasetet
        day_phasetet = {};
        for day = unique(epochs(:,1))'
            tet_numcells = zeros(99,1);
            for ep = epochs(find(epochs(:,1)==day),2)'
                
                if (~isempty(cellfilter))
                    tetlist =  evaluatefilter(cellinfo{day}{ep}, ...
                        cellfilter);
                    tetlist = unique(tetlist(:,1))';
                elseif ~isempty(tetfilter)
                    tetlist =  evaluatefilter(tetinfo{day}{ep}, ...
                        tetfilter);
                    tetlist = unique(tetlist(:,1))';
                end
                % tabulate # cells over all epochs that day
                for tet = tetlist
                    try
                        tet_numcells(tet) = tet_numcells(tet) + tetinfo{day}{ep}{tet}.numcells;
                    catch
                    end
                end
                
            end
            if any(tet_numcells)
                [dummy phasetet] = max(tet_numcells);
                day_phasetet{day} = phasetet;
            else
                day_phasetet{day} = {};
            end
        end
    end



for i = 1:size(epochs,1)
    
    
    % if onephasetet specified, just grab the one tetrode (make sure ngammal is 1 in datafilter script)
    if onephasetet
        tetlist = day_phasetet{epochs(i,1)};
    % otherwise proceed normally
    elseif (~isempty(cellfilter))
        tetlist =  evaluatefilter(cellinfo{epochs(i,1)}{epochs(i,2)}, ...
            cellfilter);
        tetlist = unique(tetlist(:,1))';
        %if tetfilter, then do that
    elseif ~isempty(tetfilter)
        tetlist =  evaluatefilter(tetinfo{epochs(i,1)}{epochs(i,2)}, ...
            tetfilter);
        tetlist = unique(tetlist(:,1))';
    end
    
    if (~isempty(tetlist))
        
        % go through the tetlist and construct an array where each element
        % represents the number of active tetrodes for each 1 ms timestep.
        
        try
            r = rip{epochs(i,1)}{epochs(i,2)}{tetlist(1)};
        catch
            keyboard
        end
        times = r.timerange(1):0.001:r.timerange(end);
        ngammal = zeros(size(times));
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
                disp(sprintf('excluded last gammal, day %d epoch %d',epochs(i,1),epochs(i,2)))
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
                ngammal = ngammal + interp1(tlist(:,1), tlist(:,2), times, 'nearest');
            catch
                keyboard
            end
        end
        out{epochs(i,1)}{epochs(i,2)}.time = times;
        clear times;
        out{epochs(i,1)}{epochs(i,2)}.ngammal = ngammal;
        
        % print out % gammal time of each epoch
        %percent_epochtime = round(100 * sum(ngammal>0)/length(ngammal));
        %    disp(sprintf('d %d e %d, %d%% gamma > 0',epochs(i,1),epochs(i,2),percent_epochtime))
        percent_epochtime = round(100 * sum(ngammal>=ngammal_report)/length(ngammal));
        disp(sprintf('d %d e %d, %d%% gamma >= %d',epochs(i,1),epochs(i,2),percent_epochtime,round(ngammal_report)))
        
        
    else    % empty tetlist this epoch
        
        disp(sprintf('getgammaltimes -- no valid tetrodes: day %d epoch %d >> ignoring epoch',epochs(i,1),epochs(i,2)))
        % retrieve some empty times structure for that epoch
        for tet=1:7
            if ~isempty(rip{epochs(i,1)}{epochs(i,2)}{tet})
                r = rip{epochs(i,1)}{epochs(i,2)}{tet};
                times = r.timerange(1):0.001:r.timerange(end);
                continue
            end
        end
        out{epochs(i,1)}{epochs(i,2)}.time = times;
        out{epochs(i,1)}{epochs(i,2)}.ngammal = zeros(size(times));
        clear times;
    end
    
end
end
