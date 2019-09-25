function [out] = gethighthetatimes2(animaldir,animalprefix, epochs, tetlist,varargin)
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
%   'powerratio1' --
        % theta:supratdelta
%   'powerratio2' --
        % theta:delta
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


% assign the options
cellfilter = '';

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

% old code
%loaddays = unique(epochs(:,1));
%rip = loaddatastruct(animaldir, animalprefix, 'ripples', loaddays);
%

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
	% represents the number of active (high theta) tetrodes for each 1 ms timestep.
	try
        theta = loadeegstruct(animaldir, animalprefix, 'theta', epochs(i,1), epochs(i,2),tetlist(1));
        delta = loadeegstruct(animaldir, animalprefix, 'delta', epochs(i,1), epochs(i,2),tetlist(1));
        supratheta = loadeegstruct(animaldir, animalprefix, 'supratheta', epochs(i,1), epochs(i,2),tetlist(1));
        pos = loaddatastruct(animaldir, animalprefix, 'pos', epochs(i,1));
        Fs = theta{epochs(i,1)}{epochs(i,2)}{tetlist(1)}.samprate;
	catch
	    keyboard
    end
    
    % establish time vector (at 150 Hz, Fs of filtered bands) from beginning to end of epoch
    times = geteegtimes(theta{epochs(i,1)}{epochs(i,2)}{tetlist(1)});
    out{epochs(i,1)}{epochs(i,2)}.time = times';     
    ntheta = zeros(size(times));
    
    for t = 1:length(tetlist)        
        % find all thresholded periods in epoch
        figure
        hold on

        %power
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
        ntheta = ntheta' + tvalid;
	end
     
    % threshold velocity for each time point (maybe quite slow, but just for testing)
    for k=1:length(ntheta)
        postimes = pos{epochs(i,1)}{epochs(i,2)}.data(:,1);
        vel = pos{epochs(i,1)}{epochs(i,2)}.data(lookup(times(k),postimes),5);
        if vel < velocitythresh
            ntheta(k)=0;
        end
    end
    
    % remove periods that don't exceed mindur
   
    % find start and end indices of thresholded periods
    startindices = find([0 diff(ntheta') 0] == 1);       % 1 of this vector 
    endindices = find([0 diff(ntheta') 0] == -1) - 1;
    indices = [startindices' endindices'];
    durations_ind = indices(:,2)-indices(:,1) + ones(size(indices,1),1);
    % cross check -- plot durations of each period
    durations_sec = durations_ind/Fs;
        figure
        hist(durations_sec,0:.1:10)
       
    threshperiods = find(durations_sec > mindur);    % row #s (in "indices" vector) of valid periods
    
	% lastly, clear ntheta and install series of 1s for each valid period
    ntheta = zeros(size(times))';               
    for j=1:length(threshperiods)
        onesvec = ones(durations_ind(threshperiods(j)),1);
        ntheta(indices(threshperiods(j),1): ...
               indices(threshperiods(j),2)) = onesvec; 
    end
    out{epochs(i,1)}{epochs(i,2)}.ntheta = ntheta;       % # of theta power ratio thresholds
    end
end


















