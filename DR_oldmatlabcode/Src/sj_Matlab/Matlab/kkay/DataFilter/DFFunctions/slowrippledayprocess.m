function slowrippledayprocess(directoryname,fileprefix,days, varargin)
%SLOWRIPPLEDAYPROCESS(directoryname,fileprefix,days, options)
%
%Applies a slow ripple filter (100-130 Hz) to all epochs for each day and 
%saves the data in in the EEG subdirectory of the directoryname folder.  
%
%directoryname - example '/data99/user/animaldatafolder', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile (e.g. 'fre')
%
%days -          a vector of experiment day numbers 
%
%options -
%		'daytetlist', [day tet ; day tet ...]
%			specifies, for each day, the tetrodes for which ripple
%			extraction should be done
%       'area_tet', Default: 'CA3'
%           specifies the area_tet for which slow ripple extraction should be
%           done. This assumes that a tetinfostruct exists.
%		'f', matfilename
%			specifies the name of the mat file containing the
%			ripple filter to use 
%			(default /home/mcarr/Src/Matlab/Filters/slowripplefilter.mat).

daytetlist = [];
f = '';
defaultfilter = '/home/mcarr/Src/Matlab/Filters/slowripplefilter.mat';

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'daytetlist'
            daytetlist = varargin{option+1};
        case 'f'
    	    f = varargin{option+1};
        case 'area'
            area_tet = varargin{option+1};
    end
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
    % create the list of files for this day that we should filter
    if isempty(daytetlist) % && isempty(area_tet) 
       tmpflist = dir(sprintf('%s/EEG/*eeg%02d-*.mat', directoryname, day));
       flist = cell(size(tmpflist));
       for i = 1:length(tmpflist)
	   flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
       end
     elseif ~isempty(area_tet)
        flist = {};
        load(sprintf('%s/%stetinfo.mat',directoryname,fileprefix));
        area_tetstring = sprintf('isequal($area_tet,''%s'')',area_tet);
        tmptetlist = evaluatefilter(tetinfo,area_tetstring);
        tet = unique(tmptetlist(tmptetlist(:,1)==day,3));
        tmpflist = [];
        for t = 1:length(tet);
            tmp = dir(sprintf('%s/EEG/*eeg%02d-*-%02d.mat', ...
                directoryname, day,tet(t)));
            tmpflist = [tmpflist; tmp];
        end
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
    temp = filtfilt(slowripplefilter,1,eeg{day}{epoch}{tet}.data');
    hdata = hilbert(temp);
    env = abs(hdata);
    phase = angle(hdata);

    slowripple{day}{epoch}{tet}.samprate = eeg{day}{epoch}{tet}.samprate;
    slowripple{day}{epoch}{tet}.starttime = eeg{day}{epoch}{tet}.starttime;

    slowripple{day}{epoch}{tet}.data(:,1) = int16(temp);
    slowripple{day}{epoch}{tet}.data(:,2) = int16(phase*10000);
    slowripple{day}{epoch}{tet}.data(:,3) = int16(env);
    slowripple{day}{epoch}{tet}.fields = ...
        'filtered_amplitude instantaneous_phase*10000 envelope_magnitude';

    clear eegrec
  	% replace the filtered invalid entries with the minimum int16 value of
	% -32768
	for i = 1:length(lo)
	    if lo(i) > 1 && hi(i) < length(slowripple{day}{epoch}{tet}.data)
            slowripple{day}{epoch}{tet}.data(lo(i):hi(i)) = minint;
	    end
	end

	% save the resulting file
	ripplefile = sprintf('%s/EEG/%sslowripple%02d-%d-%02d.mat', ...
			    directoryname, fileprefix, day, epoch, tet);
        save(ripplefile, 'slowripple');
	clear slowripple
    end
end
