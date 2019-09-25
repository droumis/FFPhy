function processdays_new(folderprefix, animdirect,fileprefix, days, varargin)
%PROCESSDAYS(FOLDERPREFIX, ANIMDIRECT, FILEPREFIX, DAYS, options)
%Run this program from the folder containing all of the animal's day folders
%It assumes that the folders are called animalname01,animalname02, ... , animalname10, ...
%FOLDERPREFIX -- if the folders are called conley01, conley02, ... then FOLDERPREFIX is 'conley'
%ANIMDIRECT -- folder where all processed data will be stored.  Example: '/data13/mkarlsso/Con/'
%FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.
%DAYS -- this is a row vector containing the day numbers to process. 
%
%OPTIONS -- optional input in the form of: 'option', value, 'option', value
%           optional arguments are:
%           DIODEPOS -- 0 uses back diode for pos, 1 uses front diode for
%               pos, and values in between use the proportionate distance
%               between the two diodes. (default 0.5, ie, average of diodes)
%           CMPERPIX -- size of each pixel in centimeters (default 1)
%           POSFILT -- filter to use for pos smoothing when computing
%               velocity (default gaussian(30*0.2, 60))
%           PROCESSEEG -- 1 to process EEG, 0 to skip EEG processing (default 1)
%           SYSTEM -- 1 for old rig, 2 for nspike (default 2)
%           VARPREFIX -- the first three letters of the animals name, which is attached to all variable names.
%              I recommend putting '' for this, which will not include any name in the variables. (Default '')
%           REVERSEX -- mirror image pos info along x axis (default 0)
%           REVERSEY -- mirror image pos info along y axis (default 0)

diodepos = 0.5;
cmperpix = 1;
processeeg = 1;
system = 2;
posfilt = gaussian(30*0.2, 60);
varprefix = '';
reversex = 0;
reversey = 0;

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



for i = days
   if i < 10
      daystring = ['0',num2str(i)];
   else
      daystring = num2str(i);
   end
   daydirect = [folderprefix,daystring];
   disp(['Day ',daystring]);
   dayprocess(daydirect, animdirect,fileprefix, i, 'diodepos', diodepos, 'cmperpix', cmperpix, 'processeeg', processeeg, 'system', system, 'posfilt', posfilt, 'varprefix', varprefix, 'reversex',reversex,'reversey',reversey);
end
