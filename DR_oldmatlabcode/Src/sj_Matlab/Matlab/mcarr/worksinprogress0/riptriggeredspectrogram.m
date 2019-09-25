function riptriggeredspectrogram(directoryname,fileprefix,days,varargin)

% out = riptriggeredspectrogram(directoryname,fileprefix,days,varargin)
%  Computes and saves the spectrogram around each ripple.
%
%directoryname - example '/data99/user/animaldatafolder', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile (e.g. 'fre')
%
%days -          a vector of experiment day numbers 
%

%  Options:
%       fpass-- Determines the frequency range for computing spectrum.
%           Default: [2 350]
%       average_trials-- Determines if events are averaged or not.
%           Default: 0
%       spectrum_window-- Determines the sliding window used to compute
%           the event triggered spectrogram. Default: [0.1 0.01]
%       event_window--Determines the size of the window around each
%           triggering event. Default: [0.2 0.4]
%       nonreference-- Specifies whether to use EEG or EEGnonreference to
%           complete filtering. Default: 1
%       cellfilter--Determines which tetrodes are used to detect triggering
%           events. Default: 'isequal($area,''CA1'')'
%       tetfilter--Determines which tetrodes ripple triggered spectrum is
%           computed for. Default: 'isequal($area,''CA1'')|isequal($area,''CA3'')'

%parse the options
params = {};
params.Fs = 1500;
params.fpass = [2 350];
params.trialave = 0;
win = [0.2 0.4];
cwin = [0.1 0.01];
daytetlist = [];
cellfilter = 'isequal($area,''CA1'')';
tetfilter = 'isequal($area,''CA1'')|isequal($area,''CA3'')';
nonreference = 1;
savedirectoryname = directoryname;

for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})       
        switch(varargin{option})
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'aveargetrials'
                params.trialave = varargin{option+1};
            case 'spectrum_window'
                cwin = varargin{option+1};
            case 'event_window'
                win = varargin{option+1};
            case 'cellfilter'
                cellfilter = varargin{option+1};
            case 'daytetlist'
                daytetlist = varargin{option+1};
            case 'nonreference'
                nonreference = varargin{option+1};
            case 'tetfilter'
                tetfilter = varargin{option+1};
            case 'savedirectory'
                savedirectoryname = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end


days = days(:)';
for day = days
    %Load ripple structure & cellinfo structures for day
    ripples = loaddatastruct(directoryname,fileprefix,'ripples',day);
    cellinfo = loaddatastruct(directoryname,fileprefix,'tetinfo',day);
    valid = evaluatefilter(cellinfo{day},tetfilter);
    valid = unique(valid(:,2));
    %Go through each run epoch and compute the ripple triggered spectrogram
    %for each valid tetrod
    for epoch = 1:min(length(ripples{day}),7)
        if ~isempty(ripples{day}{epoch})
            %Find valid riptimes
            riptimes = getripples([day epoch], ripples, cellinfo, 'cellfilter', cellfilter,'minstd',3);

            %Create list of tetrodes for this day
            if isempty(daytetlist)
                if nonreference
                    tmpflist = dir(sprintf('%s/EEGnonreference/*eeg%02d-%d-*.mat', directoryname, day,epoch));
                else
                    tmpflist = dir(sprintf('%s/EEG/*eeg%02d-%d-*.mat',directoryname,day,epoch));
                end
                flist = cell(size(tmpflist));
                for i = 1:length(tmpflist)
                    if nonreference
                        flist{i} = sprintf('%s/EEGnonreference/%s', directoryname, tmpflist(i).name);
                    else
                        flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
                    end
                end
            else
                flist = {};
                tet = daytetlist(find(daytetlist(:,1)==day),2);
                for t = tet
                    if nonreference
                        tmpflist = dir(sprintf('%s/EEGnonreference/*eeg%02d-%d-%02d.mat', ...
                            directoryname, day, epoch, t));
                    else
                        tmpflist = dir(sprintf('%s/EEGnonreference/*eeg%02d-%d-%02d.mat', ...
                            directoryname, day, epoch, t));
                    end
                    for i = 1:length(tmpflist)
                        if nonreference
                            flist{i} = sprintf('%s/EEGnonreference/%s', directoryname, tmpflist(i).name);
                        else
                            flist{i} = sprintf('%s/EEG/%s', directoryname, tmpflist(i).name);
                        end
                    end
                end
            end

            %Go through each file in flist and compute the riptriggered
            %spectrogram if that tetrode is a valid tetrode
            for fnum = 1:length(flist)
                %get the terode number
                dash = find(flist{fnum} == '-');
                tet = str2double(flist{fnum}((dash(2)+1):dash(2)+3));

                %Determine if tetrode is valid
                validtet = any(tet==valid);
                if validtet
                    load(flist{fnum});
                    e = eeg{day}{epoch}{tet}.data';
                    starttime = eeg{day}{epoch}{tet}.starttime;
                    endtime = (length(e)-1) * (1 / params.Fs);

                    % Define triggering events as the start of each ripple
                    triggers = riptimes(:,1)-starttime;

                    %Remove triggering events that are too close to the beginning or end
                    while triggers(1)<win(1)
                        triggers(1) = [];
                    end
                    while triggers(end)> endtime-win(2)
                        triggers(end) = [];
                    end

                    % Calculate the event triggered spectrogram
                    [S,t,f] = mtspecgramtrigc(e,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);

                    % Compute a z-scored spectrogram using the mean and std for the entire session
                    P = mtspecgramc(e,[cwin(1) cwin(1)],params);
                    meanP = mean(P);
                    stdP = std(P);

                    for i = 1:size(S,1)
                        for j = 1:size(S,3)
                            S(i,:,j) = (S(i,:,j) - meanP)./stdP;
                        end
                    end
                    spectrum{day}{epoch}{tet}.spectrum = S;
                    spectrum{day}{epoch}{tet}.time = t-win(1);
                    spectrum{day}{epoch}{tet}.frequency = f;
                    spectrum{day}{epoch}{tet}.ripples = triggers + starttime;

                    %Save the resulting file
                    if nonreference
                        spectrumfile = sprintf('%s/EEGnonreference/%sspectrum%02d-%d-%02d.mat', ...
                            savedirectoryname,fileprefix,day,epoch,tet);
                    else
                        spectrumfile = sprintf('%s/EEG/%sspectrum%02d-%d-%02d.mat', ...
                            savedirectoryname,fileprefix,day,epoch,tet);
                    end

                    if exist('spectrum','var')
                        save(spectrumfile,'spectrum');
                        clear spectrum
                    end
                end
            end
        end
    end
end