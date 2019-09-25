function currdir = sj_makeposfromrawpos1(rawpos, animdirect, fileprefix, daynum, epochs, varargin)


% Shantanu 07 Jan 2009
% Using dayprocess.m from /home/loren/Src/Matlab
% Adding multiunit and dio processin from Lorens and Calebs code in
% /home/loren/Src/Matlab/diomatlab

%%

%Run this program from the folder containing all of the animal's day folders.  It will process the matclust files for the desired
%day, and put the files in the ANIMDIRECT directory.
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
%           DIODEPOS -- 0 uses back diode for pos, 1 uses front diode for
%               pos, and values in between use the proportionate distance
%               between the two diodes. (default 0.5, ie, average of diodes)
%           CMPERPIX -- size of each pixel in centimeters (must be specified)
%           POSFILT -- filter to use for pos smoothing when computing
%               velocity (default gaussian(30*0.5, 60))
%           PROCESSEEG -- 1 to process EEG, 0 to skip EEG processing (default 1)
%           SYSTEM -- 1 for old rig, 2 for nspike (default 2)
%           VARPREFIX -- the first three letters of the animals name, which is attached to all variable names.
%              I recommend putting '' for this, which will not include any name in the variables. (Default '')
%           REVERSEX -- mirror image pos info along x axis (default 0)
%           REVERSEY -- mirror image pos info along y axis (default 0)

%%

% set default values for the optional arguments
diodepos = 0.5;
cmperpix = NaN;
processeeg = 1;
system = 2;
posfilt = gaussian(30*0.5, 60);
varprefix = '';
reversex = 0;
reversey = 0;
processdio = 1;


% replace default values with any supplied ones
for option = 1:2:length(varargin)-1
    
    switch varargin{option}
        case 'diodepos'
            diodepos = varargin{option+1};
        case 'cmperpix'
            cmperpix = varargin{option+1};
        case 'processeeg'
            processeeg = varargin{option+1};
        case 'system'
            system = varargin{option+1};
        case 'posfilt'
            posfilt = varargin{option+1};
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

if (~isfinite(cmperpix))
    error('cmperpix must be specified');
end

%%
currdir = pwd;
if (animdirect(end) == '/')
    animdirect = animdirect(1:end-1);
end
%cd(animdirect)
%cd(currdir)

NUMTETRODES = 30;
currdir = pwd;

if (daynum < 10)
    daystring = ['0',num2str(daynum)];
else
    daystring = num2str(daynum);
end

for epoch=epochs
    
    % interpolate over the raw positions to get location and direction
    % at each time
    if ~isempty(rawpos{daynum}{epoch}.data)
        pos{daynum}{epoch} = posinterp(rawpos{daynum}{epoch}, 'diode', diodepos,...
            'maxv', 300, 'maxdevpoints', 30, 'reversex', reversex, 'reversey', reversey);
        
        % multiply pos data by cmperpix to convert from pixel units to cm
        pos{daynum}{epoch}.cmperpixel = cmperpix;
        pos{daynum}{epoch}.data(:,2:3) = pos{daynum}{epoch}.data(:,2:3)*cmperpix;
        
        % now run functions to add additional behavior parameters to pos struct
        pos{daynum}{epoch} = addvelocity(pos{daynum}{epoch}, posfilt);
    else
        pos{daynum}{epoch}.data = [];
        pos{daynum}{epoch}.cmperpixel = cmperpix;
    end
    
end

cd(animdirect);
eval([varprefix,'pos = pos;']);
eval(['save ',fileprefix,'pos',daystring,' ',varprefix,'pos']);
cd(currdir);
%%

%if (processdio)
%    creatediostruct(daydirect,animdirect,fileprefix,daynum);
%end
%%


