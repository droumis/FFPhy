function [out] = gethighthetatimes2(animaldir,animalprefix, epochs, tetlist, varargin)
% out = gethighthetatimes(animaldir,animalprefix,epochs, tetlist, options)
%
%  Re-tooled version of getriptimes -- relies on filtered eeg bands
%  adjacent to theta, namely delta (1-4) and supratheta (12-14) a la
%  Mizuseki--Buzsaki-2009;
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
%   'powerratio' --
        % find this manually, but specified
%   'mindur' ---
        % minimum permissible high theta period  -- in SECONDS
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

smoothing_width=.2;
velocitythresh = 0;

% assign the options
cellfilter = '';
tetfilter = '';

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'statictet'
            statictet = varargin{option+1};
        case 'powerratio1'
            powerratio1 = varargin{option+1};
        case 'powerratio2'
            powerratio2 = varargin{option+1};
        case 'velocitythresh'           %% for the purposes of testing -- for actual analysis, use getlinstate or get2dstate!
            velocitythresh = varargin{option+1};
        case 'mindur'
            mindur = varargin{option+1};
        case 'tetfilter'
            tetfilter = varargin{option+1};
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

%check to see if a cell filter, or tetrode filter, is specified
if (~isempty(cellfilter))
    % this will cause us to ignore tetlist
    cellinfo = loaddatastruct(animaldir, animalprefix, 'cellinfo');
elseif (~isempty(tetfilter))
    % this will cause us to ignore tetlist
    tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo');
end

for i = 1:size(epochs,1)

    % if cellfilter is set, apply it to this day and epoch
    if ~isempty(cellfilter)
	tetlist =  evaluatefilter(cellinfo{epochs(i,1)}{epochs(i,2)}, ...
				cellfilter); 
    	% get rid of the cell indeces and extract only the tetrode numbers 
        tetlist = unique(tetlist(:,1))';
    % if tetfilter, retrieve valid tetrodes
    elseif ~isempty(tetfilter)  
	tetlist =  evaluatefilter(tetinfo{epochs(i,1)}{epochs(i,2)}, ...
				tetfilter); 
    end
    
  
	% go through the tetlist and construct an an array where each element 
	% represents the number of active (high theta) tetrodes for each
	% timestep (see below how long)
    
    % first establish time vector (at 150 Hz, Fs of filtered bands) from beginning to end of epoch

    if ~isempty(tetlist)       % only load theta if there are valid tetrodes this epoch
        theta = loadeegstruct(animaldir, animalprefix, 'theta', epochs(i,1), epochs(i,2),tetlist(1));
        times = geteegtimes(theta{epochs(i,1)}{epochs(i,2)}{tetlist(1)});
        out{epochs(i,1)}{epochs(i,2)}.time = times';     
        ntheta = zeros(size(times));
            ntheta = ntheta';
    else                       % if not, end the function with no valid theta times -- just get a times vector and a ntheta with 0s
        for tet=1:30
           try
               theta = loadeegstruct(animaldir, animalprefix, 'theta', epochs(i,1),epochs(i,2),tet);
               times = geteegtimes(theta{epochs(i,1)}{epochs(i,2)}{tet});
               if ~isempty(times)
                    out{epochs(i,1)}{epochs(i,2)}.time = times'; 
                    continue
               end
           catch
           end
        end
        disp(sprintf('no valid tetrodes for high theta times: day %d epoch %d',epochs(i,1),epochs(i,2)));
        out{epochs(i,1)}{epochs(i,2)}.ntheta = zeros(size(times))'; 
        continue 
    end
    
    pos = loaddatastruct(animaldir, animalprefix, 'pos', epochs(i,1));
    
    for t = 1:length(tetlist)       % iterate through valid tetrodes
        
        %load data
        try
            theta = loadeegstruct(animaldir, animalprefix, 'theta', epochs(i,1), epochs(i,2),tetlist(t));
                times2 = geteegtimes(theta{epochs(i,1)}{epochs(i,2)}{tetlist(t)});
            delta = loadeegstruct(animaldir, animalprefix, 'delta', epochs(i,1), epochs(i,2),tetlist(t));
            supratheta = loadeegstruct(animaldir, animalprefix, 'supratheta', epochs(i,1), epochs(i,2),tetlist(t));
            Fs = theta{epochs(i,1)}{epochs(i,2)}{tetlist(t)}.samprate;
        catch
            keyboard
        end
        
        % find all power-thresholded periods in epoch
        % powerratio1 -- theta:supratheta
        ratiotmp1 = double(theta{epochs(i,1)}{epochs(i,2)}{tetlist(t)}.data(:,3))./double(supratheta{epochs(i,1)}{epochs(i,2)}{tetlist(t)}.data(:,3));
        % powerratio2 -- theta:delta        
        ratiotmp2 = double(theta{epochs(i,1)}{epochs(i,2)}{tetlist(t)}.data(:,3))./double(delta{epochs(i,1)}{epochs(i,2)}{tetlist(t)}.data(:,3));
        % smooth the power ratios
        samprate = 150;
        kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
        ratiotmp1 = smoothvect(ratiotmp1, kernel);
        ratiotmp2 = smoothvect(ratiotmp2, kernel);
        % threshold
        tvalid = (ratiotmp1 > powerratio1) & (ratiotmp2 > powerratio2);

        % now remove periods that don't exceed mindur
        % find start and end indices of thresholded periods
        startindices = find(diff([0 tvalid']) == 1);       % 1 of this vector 
        endindices = find(diff([tvalid' 0]) == -1);
        indices = [startindices' endindices'];
        durations_ind = indices(:,2)-indices(:,1) + ones(size(indices,1),1);
        durations_sec = durations_ind/Fs;
        % cross check -- plot durations of each period
        % figure
        % hist(durations_sec,0:.1:10)
       
        threshperiods = find(durations_sec > mindur);    % row #s (in "indices" vector) of valid periods
    
        % lastly, clear tvalid and install series of ones for each valid period
        tvalid = zeros(size(times));   % adopt the length of the original times vector               
        for j=1:length(threshperiods)
            begintime = times2(indices(threshperiods(j),1));   % clocktime of beginning of specific 1s period
            endtime = times2(indices(threshperiods(j),2));      
            tvalid(lookup(begintime,times): ...
                   lookup(endtime,times)) = 1; 
        end
        
        % check possible odd error where ntheta slightly different size
        % than tvalid
               % implies that filtered eeg data is of slightly different
               % sample #s between tetrodes
        samprate=theta{epochs(i,1)}{epochs(i,2)}{tetlist(t)}.samprate;
            
        
        ntheta = ntheta + tvalid';

       
    end
    
    
    
    % threshold velocity for each time point (maybe quite slow, but just for testing)
    if velocitythresh ~= 0
        for k=1:length(ntheta)
            postimes = pos{epochs(i,1)}{epochs(i,2)}.data(:,1);
            vel = pos{epochs(i,1)}{epochs(i,2)}.data(lookup(times(k),postimes),5);
            if vel < velocitythresh
                ntheta(k)=0;
            end
        end
    end
    
    out{epochs(i,1)}{epochs(i,2)}.ntheta = ntheta;       % # of theta power ratio thresholds
    
    end
end


















