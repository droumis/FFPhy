function toppeak = eegprocess(daydirect, animdirect, fileprefix, daynum, varargin)

%Run this program from the folder containing all of the animal's day
%folders.  It only processes the eegdata and leaves all other information
%intact.
%
%DAYPROCESS(DAYDIRECT, ANIMDIRECT, VARPREFIX, FILEPREFIX, DAYNUM, OPTIONS)
%
%DAYDIRECT -- folder name where the day's data is stored
%ANIMDIRECT -- the path to where the animal's processed data will be stored -- example '/data99/student/Mil'
%FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.
%DAYNUM -- the day number for the experiment (starting with 1)
%OPTIONS -- optional input in the form of: 'option', value, 'option', value
%           optional arguments are:
%           SYSTEM -- 1 for old rig, 2 for nspike (default 2)
%           VARPREFIX -- the first three letters of the animals name, which is attached to all variable names.
%              I recommend putting '' for this, which will not include any name in the variables. (Default '')

% set default values for the optional arguments

system = 2;
varprefix = '';

% replace default values with any supplied ones
for option = 1:2:length(varargin)-1
    
    switch varargin{option}       
        case 'system'
            system = varargin{option+1};
        case 'varprefix'
            varprefix = varargin{option+1};    
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

currdir = pwd;
if (animdirect(end) == '/')
   animdirect = animdirect(1:end-1);
end
cd(animdirect)
eegdir = dir('EEG');
if (isempty(eegdir))
   %an EEG folder needs to be created
   !mkdir EEG
end
cd(currdir)

NUMTETRODES = 30;
currdir = pwd;
cd(daydirect);
times = gettimes('times');

if (daynum < 10)
   daystring = ['0',num2str(daynum)];
else
   daystring = num2str(daynum);
end

% get a list of all the eeg files in the directory
eegfiles = dir('*.eeg');
    
if ~isempty(eegfiles) % if there are eeg files to process, then proceed
    for eegfile = 1:length(eegfiles)
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
        % now get the eeg data for each epoch
        for epoch = 1:length(times)
            if ~isempty(times(epoch).starttime)
            	['epoch ', num2str(epoch),' tetrode ',num2str(eegtet)];
                starttime = times(epoch).starttime;
                endtime = times(epoch).endtime;
                if (system == 1)
                	eegstruct = readeeg(eegfilename, starttime, endtime, eegtet/NUMTETRODES);
                elseif (system == 2)
                	eegstruct = nreadeeg(eegfilename, starttime, endtime);
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