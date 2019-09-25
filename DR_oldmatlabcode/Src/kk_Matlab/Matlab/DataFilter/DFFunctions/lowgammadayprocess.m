function kk_lowgammadayprocess(directoryname,fileprefix,days, varargin)
%LOWGAMMADAYPROCESS(directoryname,fileprefix,days, options)
%
%Applies a low gamma filter to all epochs for each day and saves the data
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
%
%		'daytetlist', [day tet ; day tet ...]
%			specifies, for each day, the tetrodes for which low gamma
%			extraction should be done
%       'tetfilter', 'isequal($area,''CA1'')'
%           specifies the filter to use to determine which tetrodes 
%           low gamma extraction should be done. This assumes that a
%           tetinfostruct exists.
%		'f', filter
%			specifies the filter to use. This should be made specificially
%			for each animal based on individual cutoffs for low and high 
%           gamma.
%		'assignphase', 0 or 1
%			specifices whether to ignore spike fields (0) or assign 
%			a low gamma phase to each spike and save the data (1).
%			Default 0
%       'nonreference', 0 or 1
%           specifies whether to use EEG or EEGnonreferencehighgammasextract.m
%           Default 0

daytetlist = [];
f = '';
defaultfilter = (['/home/kkay/matlab/Filters/bonlowgammafilter.mat']);
assignphase = 0;
tetfilter = '';
reftognd = 0;
savedirectoryname = directoryname;

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'daytetlist'
            daytetlist = varargin{option+1};
        case 'f'
    	    f = varargin{option+1};
        case 'assignphase'
            assignphase = varargin{option+1};
        case 'tetfilter'
            tetfilter = varargin{option+1};
        case 'reftognd'
            reftognd = varargin{option+1};            
        case 'nonreference'
            subtractreference = varargin{option+1};
        case 'savedirectory'
            savedirectoryname = varargin{option+1};
    end
end


% check to see if the directory has a trailing '/'
if (directoryname(end) ~= '/')
    warning('directoryname should end with a ''/'', appending one and continuing');
    directoryname(end+1) = '/';
end

minint = -32768;
days = days(:)';

% if the filter was not specified, load the default
if isempty(f)
    load(defaultfilter);
else
    load(f)
    %eval(['load',f])   %older code. commented out.
end

for day = days
    if (assignphase)
        %load up the spike file
        spikes = loaddatastruct(directoryname,fileprefix, 'spikes',day);
    end
    % create the list of files for this day that we should filter
     if isempty(daytetlist) && isempty(tetfilter)
        if reftognd
            tmpflist = dir(sprintf('%s/EEG/*eeggnd%02d-*.mat', directoryname, day));
        else
            tmpflist = dir(sprintf('%s/EEG/*eeg%02d-*.mat', directoryname, day));
        end
        flist = cell(size(tmpflist));
        for i = 1:length(tmpflist)
            if reftognd
                flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
            else
                flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
            end
        end
%         keyboard
    elseif ~isempty(tetfilter)
        flist = {};
        load(sprintf('%s/%stetinfo.mat',directoryname,fileprefix));
        tmptetlist = evaluatefilter(tetinfo,tetfilter);
        tet = unique(tmptetlist(tmptetlist(:,1)==day,3));
        tmpflist = [];
        for t = 1:length(tet);
            if reftognd
                tmp = dir(sprintf('%s/EEG/*eeggnd%02d-*-%02d.mat',...
                    directoryname,day, tet(t)));
            else
                tmp = dir(sprintf('%s/EEG/*eeg%02d-*-%02d.mat', ...
                    directoryname, day,tet(t)));
            end
            tmpflist = [tmpflist; tmp];
        end
        for i = 1:length(tmpflist)
            if reftognd
                flist{i} = sprintf('%s/EEG/%s',directoryname, tmpflist(i).name);
            else
                flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
            end
        end
    else
        % find the rows associated with this day
        flist = {};
        tmpflist = [];
        tet = daytetlist(find(daytetlist(:,1) == day),2);
        for t = 1:length(tet);
            if reftognd
            tmp = dir(sprintf('%s/EEG/*eeggnd%02d-*-%02d.mat',...
                    directoryname,day, tet(t)));
            else
                tmp = dir(sprintf('%s/EEG/*eeg%02d-*-%02d.mat', ...
                    directoryname, day,tet(t)));
            end
            tmpflist = [tmpflist; tmp];
        end
        for i = 1:length(tmpflist)
            if reftognd
                flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
            else
                flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
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
        if reftognd
            entry = eeggnd{day}{epoch}{tet};
        else
            entry = eeg{day}{epoch}{tet};
        end
        a = find(entry.data < -30000);
        [lo,hi]= findcontiguous(a);  %find contiguous NaNs
        for i = 1:length(lo)
        	if lo(i) > 1 & hi(i) < length(entry.data)
                fill = linspace(entry.data(lo(i)-1), ...
                    entry.data(hi(i)+1), hi(i)-lo(i)+1);
                entry.data(lo(i):hi(i)) = fill;
            end
        end

        % filter and save the result as int16
        temp = filtfilt(lowgammafilter,1,double(entry.data'));
                                        % added double() for Anna's data
                                        % 12.28.12
        hdata = hilbert(temp);
        env = abs(hdata);
        phase = angle(hdata);
        
        if reftognd 
        lowgammagnd{day}{epoch}{tet}.samprate = eeg{day}{epoch}{tet}.samprate;
        lowgammagnd{day}{epoch}{tet}.starttime = eeg{day}{epoch}{tet}.starttime;

        lowgammagnd{day}{epoch}{tet}.data(:,1) = int16(temp);
        lowgammagnd{day}{epoch}{tet}.data(:,2) = int16(phase*10000);
        lowgammagnd{day}{epoch}{tet}.data(:,3) = int16(env);
        lowgammagnd{day}{epoch}{tet}.fields = ...
            'filtered_amplitude instantaneous_phase*10000 envelope_magnitude';
        else
        lowgamma{day}{epoch}{tet}.samprate = eeg{day}{epoch}{tet}.samprate;
        lowgamma{day}{epoch}{tet}.starttime = eeg{day}{epoch}{tet}.starttime;

        lowgamma{day}{epoch}{tet}.data(:,1) = int16(temp);
        lowgamma{day}{epoch}{tet}.data(:,2) = int16(phase*10000);
        lowgamma{day}{epoch}{tet}.data(:,3) = int16(env);
        lowgamma{day}{epoch}{tet}.fields = ...
            'filtered_amplitude instantaneous_phase*10000 envelope_magnitude';
        end
                    
        clear eeg
        clear eeggnd
        %replace the filtered invalid entries with the minimum int16 value of -32768
        for i = 1:length(lo)
        	if lo(i) > 1 && hi(i) < length(lowgamma{day}{epoch}{tet}.data)
                if reftognd
            	lowgammagnd{day}{epoch}{tet}.data(lo(i):hi(i)) = minint;
                else
            	lowgamma{day}{epoch}{tet}.data(lo(i):hi(i)) = minint;
                end
            end
        end

        % save the resulting file
        if reftognd
            lowgammagndfile = sprintf('%sEEG/%slowgammagnd%02d-%d-%02d.mat', ...
            	savedirectoryname, fileprefix, day, epoch, tet);
            save(lowgammagndfile, 'lowgammagnd');
            clear lowgammagnd
        else
            lowgammafile = sprintf('%sEEG/%slowgamma%02d-%d-%02d.mat', ...
            	savedirectoryname, fileprefix, day, epoch, tet);
            save(lowgammafile, 'lowgamma');
            clear lowgamma
        end
    end

end
