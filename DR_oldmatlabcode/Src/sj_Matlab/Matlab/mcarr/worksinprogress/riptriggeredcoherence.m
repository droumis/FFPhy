function riptriggeredcoherence(directoryname,fileprefix,days,varargin)

% out = riptriggeredcoherence(directoryname,fileprefix,days,varargin)
%  Computes and saves the coherence around each ripple.
%
%directoryname - example '/data99/user/animaldatafolder', a folder 
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile (e.g. 'fre')
%
%days -          a vector of experiment day numbers 
%

%  Options:
%       fpass-- Determines the frequency range for computing coherence.
%           Default: [2 350]
%       average_trials-- Determines if events are averaged or not.
%           Default: 0
%       coherence_window-- Determines the sliding window used to compute
%           the event triggered coherence. Default: [0.1 0.01]
%       event_window--Determines the size of the window around each
%           triggering event. Default: [0.2 0.4]
%       nonreference-- Specifies whether to use EEG or EEGnonreference to
%           complete filtering. Default: 1
%       cellfilter--Determines which tetrodes are used to detect triggering
%           events. Default: 'isequal($area,''CA1'')'
%       tetfilter--Determines which tetrodes are used to compute coherence.
%           Default:{{'isequal($area,''CA1'')'};{'isequal($area,''CA3'')'}}

%parse the options
params = {};
params.Fs = 1500;
params.fpass = [2 350];
params.trialave = 0;
win = [0.2 0.4];
cwin = [0.1 0.01];
tetfilter = [{'isequal($area,''CA1'')'};{'isequal($area,''CA3'')'}];
cellfilter = 'isequal($area,''CA1'')';
nonreference = 1;
savedirectoryname = directoryname;

for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})       
        switch(varargin{option})
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'aveargetrials'
                params.trialave = varargin{option+1};
            case 'coherence_window'
                cwin = varargin{option+1};
            case 'event_window'
                win = varargin{option+1};
            case 'tetfilter'
                tetfilter = varargin{option+1};
            case 'cellfilter'
                cellfilter = varargin{option+1};
            case 'nonreference'
                nonreference = varargin{option+1};
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
    cellinfo = loaddatastruct(directoryname,fileprefix,'cellinfo',day);
    valid1 = evaluatefilter(cellinfo{day},tetfilter{1});
    valid1 = unique(valid1(:,2));
    valid2 = evaluatefilter(cellinfo{day},tetfilter{2});
    valid2 = unique(valid2(:,2));
    %Go through each run epoch and compute the ripple triggered spectrogram
    %for each valid tetrode
    for epoch = 1:min(7,length(ripples{day}))
        if ~isempty(ripples{day}{epoch})
            %Find valid riptimes
            riptimes = getripples([day epoch], ripples, cellinfo, 'cellfilter', cellfilter,'minstd',3);

            %Create list of tetrodes for this day
            flist1 = cell(length(valid1),1);
            for tet = 1:length(valid1)
                t = valid1(tet);
                if nonreference
                    tmpflist = dir(sprintf('%s/EEGnonreference/*eeg%02d-%d-%02d.mat', ...
                        directoryname, day, epoch, t));
                    flist1{tet} = sprintf('%s/EEGnonreference/%s', directoryname, tmpflist.name);
                else
                    tmpflist = dir(sprintf('%s/EEGnonreference/*eeg%02d-%d-%02d.mat', ...
                        directoryname, day, epoch, t));
                    flist1{tet} = sprintf('%s/EEG/%s', directoryname, tmpflist.name);
                end
            end
            flist2 = cell(length(valid2),1);
            for tet = 1:length(valid2)
                t = valid2(tet);
                if nonreference
                    tmpflist = dir(sprintf('%s/EEGnonreference/*eeg%02d-%d-%02d.mat', ...
                        directoryname, day, epoch, t));
                    flist2{tet} = sprintf('%s/EEGnonreference/%s', directoryname, tmpflist.name);
                else
                    tmpflist = dir(sprintf('%s/EEGnonreference/*eeg%02d-%d-%02d.mat', ...
                        directoryname, day, epoch, t));
                    flist2{tet} = sprintf('%s/EEG/%s', directoryname, tmpflist.name);
                end
            end

            %Go through each file in flist1 and compute the riptriggered
            %coherence of that tetrode with all of the corresponding tetrodes
            for fnum = 1:length(flist1)
                %get the terode number
                dash = find(flist1{fnum} == '-');
                tet1 = str2double(flist1{fnum}((dash(2)+1):dash(2)+3));

                if nonreference
                    coherencefile = sprintf('%s/EEGnonreference/%scoherence%02d-%d-%02d.mat', ...
                        savedirectoryname,fileprefix,day,epoch,tet1);

                else
                    coherencefile = sprintf('%s/EEG/%scoherence%02d-%d-%02d.mat', ...
                        savedirectoryname,fileprefix,day,epoch,tet1);
                end

                %Load the saved coherence file if it exists
                if exist(coherencefile,'file')
                    load(coherencefile)
                end

                load(flist1{fnum});
                e1eeg = eeg{day}{epoch}{tet1};
                e1times = geteegtimes(e1eeg);

                for fnum2 = 1:length(flist2)
                    dash = find(flist2{fnum2} == '-');
                    tet2 = str2double(flist2{fnum2}((dash(2)+1):dash(2)+3));
                    %Compute coherence only when tet1 and tet2 are different
                    if tet1 ~= tet2
                        load(flist2{fnum2})
                        e2 = eeg{day}{epoch}{tet2};
                        e2times = geteegtimes(e2);

                        if length(e1times)>length(e2times)
                            temp = lookup(e2times,e1times);
                            e1times = e1times(temp);
                            e1 = e1eeg.data(temp);
                            e2 = e2.data;
                        elseif length(e2times)>length(e1times)
                            temp = lookup(e1times,e2times);
                            e1 = e1eeg.data;
                            e2 = e2.data(temp);
                        elseif length(e1times)==length(e2times)
                            e1 = e1eeg.data;
                            e2 = e2.data;
                        end
                        starttime = e1times(1);
                        endtime = (length(e1)-1) * (1 / params.Fs);

                        % Define triggering events as the start of each ripple
                        triggers = riptimes(:,1)-starttime;

                       %Remove triggering events that are too close to the beginning or end
                        while triggers(1)<win(1) + 0.5
                            triggers(1) = [];
                        end
                        while triggers(end)> endtime-win(2)-0.5
                            triggers(end) = [];
                        end

                        % Calculate the event triggered coherence
                        data1 = createdatamatc(e1,triggers,params.Fs,[win(1) win(2)]);
                        data2 = createdatamatc(e2,triggers,params.Fs,[win(1) win(2)]);

                        [C,phi,S12,S1,S2,t,f] = cohgramc(data1,data2,[cwin(1) cwin(2)],params);

                        coherence{day}{epoch}{tet1}{tet2}.coherence = im2uint8(C);
                        coherence{day}{epoch}{tet1}{tet2}.phase = phi;
                        coherence{day}{epoch}{tet1}{tet2}.time = t-win(1);
                        coherence{day}{epoch}{tet1}{tet2}.frequency = f;
                        coherence{day}{epoch}{tet1}{tet2}.ripples = triggers+starttime;
                    end
                end
                %Save the resulting file
                if exist('coherence','var')
                    save(coherencefile,'coherence');
                    clear coherence
                end
            end
        end
    end
end