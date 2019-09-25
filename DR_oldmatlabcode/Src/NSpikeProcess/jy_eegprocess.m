function toppeak = jy_eegprocess(animdirect, daydirect, fileprefix, daynum, varargin)
%It will process the matclust files for eeg for the desired
%day, and put the files in the ANIMDIRECT directory.
%
%based on parts of DAYPROCESS(DAYDIRECT, ANIMDIRECT, VARPREFIX, FILEPREFIX, DAYNUM, OPTIONS)
%
%DAYDIRECT  -- folder with raw data files in /data14/jai/
%              no extra / needed eg. 'H2'
%ANIMDIRECT -- folder with the animals processed data in /data14/jai/
%              no extra / needed eg. 'H2_'
%FILEPREFIX -- name of animal eg. 'H2'
%DAYNUM     -- the day number for the experiment (starting with 1)
%OPTIONS    -- optional input in the form of: 'option', value, 'option', value
%
%           SYSTEM -- 1 for old rig, 2 for nspike (default 2)
%           VARPREFIX -- the first three letters of the animals name, which is attached to all variable names.
%              I recommend putting '' for this, which will not include any name in the variables. (Default '')


currdir = pwd;
% set root directory
rootdir='/data14/jai/';

% set default values for the optional arguments
processeeg = 1;
system = 2;
varprefix = '';
reversex = 0;
reversey = 0;

% replace default values with any supplied ones
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'processeeg'
            processeeg = varargin{option+1};
        case 'system'
            system = varargin{option+1};
        case 'varprefix'
            varprefix = varargin{option+1};
        case 'reversex'
            reversex = varargin{option+1};
        case 'reversey'
            reversey = varargin{option+1};
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

% specify raw and processed data directories
if (daynum < 10)
    daystring = ['0',num2str(daynum)];
else
    daystring = num2str(daynum);
end


daydirect=strcat(rootdir,daydirect,'/',fileprefix,daystring,'/');

animdirect=strcat(rootdir,animdirect,'/');

cd(animdirect)

eegdir = dir('EEG');
if (isempty(eegdir))
    %an EEG folder needs to be created
    !mkdir EEG
end

cd(currdir)

NUMTETRODES = 30;
currdir = pwd;

% get time file for day
tfilename = strcat(daydirect,'times.mat');
times = gettimes(tfilename);

runepochs = [];
eeg = [];

%
cd(daydirect);
% get a list of all the eeg files in the directory
eegfiles = dir('*.eeg');
if (isempty(eegfiles))
    % No *.eeg, so look for *.cont
    cont_flag=1;
    eegfiles = dir('*.cont');
    
    % Check if files can be opened
    contfile = eegfiles(1).name;
    fid = fopen(contfile);
    if (fid == -1)
        error(sprintf('Cont File %s cannot be opened - Set processeeg to 0 and re-run', contfile));
    end
    
    % Get Sampling Rate and nsamplesperbuffer
    headerstr = '';
    while (~strncmp(headerstr, '%%ENDHEADER', 11))
        headerstr = fgets(fid);
        ind = strfind(headerstr, 'nsamplesperbuf');
        slen = length('nsamplesperbuf');
        if (~isempty(ind))
            nsamplesperbuf = str2num(headerstr((ind+slen):end))
        end
        ind = strfind(headerstr, 'samplingrate');
        slen = length('samplingrate');
        if (~isempty(ind))
            samplingrate = str2num(headerstr((ind+slen):end))
        end
    end
    fclose(fid);
    
    if samplingrate > 2000
        error(sprintf('Sampling rate of %d is too high - Set processeeg to 0 and re-run', samplingrate));
    end
    
end

if (processeeg) %this part filters the eeg data
    
    for eegfile = 1:length(eegfiles)
        if ~isempty(eegfiles) % if there are eeg files to process, then proceed
            
            
            % get the tetrode number and depth, which is contained in the
            % name of the eeg file
            eegfilename = eegfiles(eegfile).name;
            dashfind = strfind(eegfilename,'-');
            pointfind = strfind(eegfilename,'.');
            depth = str2num(eegfilename(dashfind(1)+1:pointfind-1));
            eegtet = str2num(eegfilename(1:dashfind-1));
            % prepare a string label for this eeg data
            if (eegtet < 10)
                eegstring = ['0',num2str(eegtet)];
            else
                eegstring = num2str(eegtet);
            end
            fid = fopen(eegfilename);
            if (fid == -1)
                error(sprintf('Cont File %s cannot be opened - Set processeeg to 0 and re-run', eegfilename));
            end
            % Get nsamples per buffer and sampling rate for each eeg file separately. This is important, since the
            % values can differ for DSPs. Updated 15 Jun 2012. SJ.
            
            % Get Sampling Rate and nsamplesperbuffer
            headerstr = '';
            while (~strncmp(headerstr, '%%ENDHEADER', 11))
                headerstr = fgets(fid);
                ind = strfind(headerstr, 'nsamplesperbuf');
                slen = length('nsamplesperbuf');
                if (~isempty(ind))
                    nsamplesperbuf = str2num(headerstr((ind+slen):end));
                end
                ind = strfind(headerstr, 'samplingrate');
                slen = length('samplingrate');
                if (~isempty(ind))
                    samplingrate = str2num(headerstr((ind+slen):end));
                end
            end
            fclose(fid);
            
            
            
            
            % now get the eeg data for each epoch
            for epoch = 1:length(times)
                if ~isempty(times(epoch).starttime)
                    ['epoch ', num2str(epoch),' tetrode ',num2str(eegtet)]
                    starttime = times(epoch).starttime;
                    endtime = times(epoch).endtime;
                    epochname = times(epoch).name;
                    if (system == 1)
                        eegstruct = readeeg(eegfilename, starttime, endtime, eegtet/NUMTETRODES);
                    elseif (system == 2)
                        if cont_flag==1
                            eegstruct = nreadcont(eegfilename, samplingrate, nsamplesperbuf, starttime, endtime);
                        else
                            eegstruct = nreadeeg(eegfilename, starttime, endtime);
                        end
                    end
                    eegstruct.depth = depth;
                    eval([varprefix,'eeg{',num2str(daynum),'}{',num2str(epoch),'}{',num2str(eegtet),'} = eegstruct;']);
                    tmpfilename = [animdirect,'/EEG/',fileprefix,'eeg',daystring,'-',num2str(epoch),'-',eegstring];
                    eval(['save ',tmpfilename, ' ',varprefix,'eeg']);
                    eval(['clear ',varprefix,'eeg']);
                end
            end
        end
    end
end




