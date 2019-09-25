function makethetastruct(directoryname,fileprefix, days, varargin) %MAKETHETASTRUCT(directoryname,fileprefix,days, options)
%
% Creates three theta structures, thetawav, thetaphase, and thetaamp which
% contain the waveform, phase and amplitude of the theta for the specificied
% tetrodes.  A separate set of structures is saved for each day
%
%directoryname - example '/data99/user/animaldatafolder/', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile (e.g. 'fre')
%
%days -          a vector of experiment day numbers 
%
%options -
%		'tetlist', [tet tet ...]
%			specifies, the set of tetrodes whose data should be
%			included.  Using this option is recommended as the
%			structures can get quite large.
%

tetlist = 1:100;
%set variable options
tetlist = [];
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'tetlist'
	    tetlist = varargin{option+1};
    end
end

for d = days
   ([animaldir, aname, 'te

j
minint = -32768;

% if the filter was not specified, load the default
if isempty(f) 
    eval(['load ', defaultfilter]);
else
    eval(['load ', f]);
end

days = days(:)';

for day = days
    % create the list of files for this day that we should filter
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
	theta{day}{epoch}{tet} = filtereeg2(eeg{day}{epoch}{tet}, ...
			thetafilter, 'int16', 1); 
        % downsample if requested 
	if (downsample ~= 1)
	    theta{day}{epoch}{tet}.data = ...
		theta{day}{epoch}{tet}.data(1:downsample:end, :);
	    theta{day}{epoch}{tet}.samprate = ...
	    	theta{day}{epoch}{tet}.samprate / downsample;
	end
	clear eegrec
	% replace the filtered invalid entries with the minimum int16 value of
	% -32768
	for i = 1:length(lo)
	    if lo(i) > 1 & hi(i) < length(theta{day}{epoch}{tet}.data)
		theta{day}{epoch}{tet}.data(lo(i):hi(i)) = minint;
	    end
	end

	% save the resulting file
	thetafile = sprintf('%s/EEG/%stheta%02d-%d-%02d.mat', ...
			    directoryname, fileprefix, day, epoch, tet);
        save(thetafile, 'theta');
	clear theta
    end
end
