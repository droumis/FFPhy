function flipeegdays(directoryname,fileprefix, days, varargin)
%FLIPEEGDAYS(ANIMDIRECT, FILEPREFIX, DAYS, options)

%Multiplies all eeg traces by -1 and saves the data in EEG and if in the
%EEGnonreference folder if applicable
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
%       'nonreference', 0 or 1
%           specifies whether to check for EEGnonreference files
%           Default 0

nonreference = 0;
daytetlist = [];

% replace default values with any supplied ones
for option = 1:2:length(varargin)-1
    
    switch varargin{option}
        case 'nonreference'
            nonreference = varargin{option+1};
        case 'daytetlist'
            daytetlist = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end


days = days(:)';

for day = days
    % create the list of files for this day that we should load
    if (isempty(daytetlist))
        tmpflist = dir(sprintf('%s/EEG/*eeg%02d-*.mat', directoryname, day));
        flist = cell(size(tmpflist));
        for i = 1:length(tmpflist)
        	flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
        end
    else
    	% find the rows associated with this day
        flist = {};
        tet = daytetlist(find(daytetlist(:,1) == day),2);
        for t = tet;
            tmpflist = dir(sprintf('%s/EEG/*eeg%02d-*-%02d.mat', ...
            	directoryname, day, t));
            for i = 1:length(tmpflist)
               flist{iind} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
            end
        end
    end
    
    % go through each file in flist and multiply by -1
    for fnum = 1:length(flist)
        % get the tetrode number and epoch
        dash = find(flist{fnum} == '-');
        epoch = str2num(flist{fnum}((dash(1)+1):(dash(2)-1)));
        tet = str2num(flist{fnum}((dash(2)+1):(dash(2)+3)));

        %load the eeg file
        load(flist{fnum});
        eeg{day}{epoch}{tet}.data = -1*eeg{day}{epoch}{tet}.data;
    
        eegfile = sprintf('%s/EEG/%seeg%02d-%d-%02d.mat',...
            directoryname,fileprefix,day,epoch,tet);
        save(eegfile,'eeg');
    end
    
    % If EEGnonreference is specified, multiply eeg files in
    % EEGnonreference by -1 as well.
    if nonreference
        % create the list of files for this day that we should load
        if (isempty(daytetlist))
            tmpflist = dir(sprintf('%s/EEGnonreference/*eeg%02d-*.mat', directoryname, day));
            flist = cell(size(tmpflist));
            for i = 1:length(tmpflist)
                flist{i} = sprintf('%s/EEGnonreference/%s', directoryname, tmpflist(i).name);
            end
        else
            % find the rows associated with this day
            flist = {};
            tet = daytetlist(find(daytetlist(:,1) == day),2);
            for t = tet;
                tmpflist = dir(sprintf('%s/EEGnonreference/*eeg%02d-*-%02d.mat', ...
                    directoryname, day, t));
                for i = 1:length(tmpflist)
                   flist{iind} = sprintf('%s/EEGnonreference/%s', directoryname, tmpflist(i).name);
                end
            end
        end

        % go through each file in flist and multiply by -1
        for fnum = 1:length(flist)
            % get the tetrode number and epoch
            dash = find(flist{fnum} == '-');
            epoch = str2num(flist{fnum}((dash(1)+1):(dash(2)-1)));
            tet = str2num(flist{fnum}((dash(2)+1):(dash(2)+3)));

            %load the eeg file
            load(flist{fnum});
            eeg{day}{epoch}{tet}.data = -1*eeg{day}{epoch}{tet}.data;

            eegfile = sprintf('%s/EEGnonreference/%seeg%02d-%d-%02d.mat',...
                directoryname,fileprefix,day,epoch,tet);
            save(eegfile,'eeg');
        end
    end
end