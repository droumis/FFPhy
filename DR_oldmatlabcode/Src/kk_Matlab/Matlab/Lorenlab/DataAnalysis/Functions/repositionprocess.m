function repositionprocess(animdirect, fileprefix, epochs, varargin)

%repositionprocess(animdirect, fileprefix, epochs, varargin)
%
%ANIMDIRECT -- the path to where the animal's processed data is stored -- example '/data99/student/Mil'
%FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.
%EPOCHS -- a list of day/epochs (Nx2) to be reprocessed
%OPTIONS -- optional input in the form of: 'option', value, 'option', value
%           optional arguments are:
%           DIODEPOS -- 0 uses back diode for pos, 1 uses front diode for
%               pos, and values in between use the proportionate distance
%               between the two diodes. (default 0.5, ie, average of diodes)
%           CMPERPIX -- size of each pixel in centimeters (default 1)
%           POSFILT -- filter to use for pos smoothing when computing
%               velocity (default gaussian(30*0.5, 60))
%           VARPREFIX -- the first three letters of the animals name, which is attached to all variable names.
%              I recommend putting '' for this, which will not include any name in the variables. (Default '')
%           REVERSEX -- mirror image pos info along x axis (default 0)
%           REVERSEY -- mirror image pos info along y axis (default 0)

% set default values for the optional arguments
diodepos = 0.5;
cmperpix = 1;
posfilt = gaussian(30*0.5, 60);
varprefix = '';
reversex = 0;
reversey = 0;

[otherOptions] = procOptions(varargin);

days = unique(epochs(:,1));

rawpos = loaddatastruct(animdirect, fileprefix, 'rawpos', days);

% load existing pos struct - will be overwritten in proper places
pos = loaddatastruct(animdirect, fileprefix, 'pos', days);

for e = 1:length(epochs)
    daynum = epochs(e,1);
    epoch = epochs(e,2);
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

savedatastruct(pos,animdirect, fileprefix, 'pos');
