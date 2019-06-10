function gammadayprocess(directoryname,fileprefix,days, varargin)
%GAMMADAYPROCESS(directoryname,fileprefix,days, options)
%
%Applies a gamma filter to all epochs for each day and saves the data in
%in the EEG subdirectory of the directoryname folder.  
%
%directoryname - example '/data99/user/animaldatafolder/', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile (e.g. 'fre')
%
%days -          a vector of experiment day numbers 
%
%options -
%		'system', 1 or 2
%			specifies old (1) or new/nspike (2) rigs. Default 2.
%
%		'daytetlist', [day tet ; day tet ...]
%			specifies, for each day, the tetrodes for which gamma
%			extraction should be done and for which gamma phase
%			should be assigned to spikes
%
%		'f', matfilename
%			specifies the name of the mat file containing the
%			gamma filter to use 
%			(default /usr/local/filtering/gammafilter.mat).  
%			Note that the filter must be called 'gammafilter'.
%		'filtergamma', 0 or 1
%			specifices whether to extract gamma from the eeg files
%			(1) or to use preexisting gamma eeg files (0).
%			Default 1;
%		'assignphase', 0 or 1
%			specifices whether to ignore spike fiels (0) or assign 
%			a gamma phase to each spike and save the data (1).
%			Default 0

daytetlist = [];
f = '';
defaultfilter = 'gammadayprocess_filter.mat';
system = 2;
filtergamma = 1;
assignphase = 0;

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'daytetlist'
	    daytetlist = varargin{option+1};
        case 'f'
	    f = varargin{option+1};
        case 'filtergamma'
	    filtergamma = varargin{option+1};
        case 'assignphase'
	    assignphase = varargin{option+1};
    end
end

% check to see if the directory has a trailing '/'
if (directoryname(end) ~= '/')
    warning('directoryname should end with a ''/'', appending one and continuing');
    directoryname(end+1) = '/';
end

minint = -32768;

% if the filter was not specified, load the default
if isempty(f) 
    eval(['load ', defaultfilter]);
else
    eval(['load ', f]);
end

days = days(:)';

for day = days
    if (assignphase)
	% load up the spike file if it exists
	spikes = loaddatastruct(directoryname, fileprefix, 'spikes', day);
    end
    % create the list of files for this day that we would filter if required
    if (isempty(daytetlist)) 
       tmpflist = dir(sprintf('%s/EEG/*eeg%02d-*.mat', directoryname, day));
       flist = cell(size(tmpflist));
       for i = 1:length(tmpflist)
	   flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
       end
    else
       % find the rows associated with this day
       flist = {};
       ind = 1;
       tet = daytetlist(find(daytetlist(:,1) == day),2);
       for t = tet;
	   tmpflist = dir(sprintf('%s/EEG/*eeg%02d-*-%02d.mat', ...
			  directoryname, day, t));
	   nfiles = length(tmpflist); 
	   for i = 1:length(tmpflist)
	       flist{iind} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
	   end
       end
    end
    
    % go through each file in flist and filter it
    for fnum = 1:length(flist)
	% get the tetrode number and epoch
	% this is ugly, but it works
	dash = find(flist{fnum} == '-');
	epoch = str2num(flist{fnum}((dash(1)+1):(dash(2)-1)));
	tet = str2num(flist{fnum}((dash(2)+1):(dash(2)+3)));
	if (filtergamma)
	    %load the eeg file
	    load(flist{fnum});
	    a = find(eeg{day}{epoch}{tet}.data < -30000);
	    [lo,hi]= findcontiguous(a);  %find contiguous NaNs
	    for i = 1:length(lo)
		if lo(i) > 1 & hi(i) < length(eeg{day}{epoch}{tet}.data)
		    fill = linspace(eeg{day}{epoch}{tet}.data(lo(i)-1), ...
			    eeg{day}{epoch}{tet}.data(hi(i)+1), hi(i)-lo(i)+1);
		    eeg{day}{epoch}{tet}.data(lo(i):hi(i)) = fill;
		end
	    end
	    % filter it and save the result as int16
	    gamma{day}{epoch}{tet} = filtereeg2(eeg{day}{epoch}{tet}, ...
			    gammafilter, 'int16', 1); 
	    clear eegrec
	    % replace the filtered invalid entries with the minimum int16 value of
	    % -32768
	    for i = 1:length(lo)
		if lo(i) > 1 & hi(i) < length(gamma{day}{epoch}{tet}.data)
		    gamma{day}{epoch}{tet}.data(lo(i):hi(i)) = minint;
		end
	    end
	    % save the resulting file
	    gammafile = sprintf('%s/EEG/%sgamma%02d-%d-%02d.mat', ...
				directoryname, fileprefix, day, epoch, tet);
	    save(gammafile, 'gamma');
	else
	    % load up the gamma file
	    gammafile = sprintf('%s/EEG/%sgamma%02d-%d-%02d.mat', ...
				directoryname, fileprefix, day, epoch, tet);
	    load(gammafile);
	end
	if (assignphase && ~isempty(spikes))
	    % check to see if there are spikes on this tetrode
	    s = [];
	    try 
		s = spikes{day}{epoch}{tet};
	    catch
	    end
	    if (~isempty(s))
		% make a list of time indeces for the gamma phases
		g = gamma{day}{epoch}{tet};
		gtimes = g.starttime:(1/g.samprate):(g.starttime + ...
					    (length(g.data)-1) / g.samprate);
		for c = 1:length(s)
		    data = [];
		    try 
			data = s{c}.data;
		    catch
		    end
		    if (length(data))
			% look up the gamma phase for each spike
			ind = lookup(data(:,1), gtimes);
			% note that the phases are * 10000;
			spikes{day}{epoch}{tet}{c}.gammaphase = ...
				double(g.data(ind,2)) / 10000;
		    end
		end
	    end
	end
	clear gamma
    end
    if (assignphase)
	% save the spike file 
	savedatastruct(spikes, directoryname, fileprefix, 'spikes');
    end
end
