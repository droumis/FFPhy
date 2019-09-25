function addtaskinfo(animdirect, fileprefix, days, epochs, varargin)

% ADDTASKINFO(animdirect, fileprefix, days, epochs, property1, value1, property2, value2, ...)
%
% This function adds descriptive information to the 'task' structure of a
% given animal. It saves the information for you.
%
% ANIMDIRECT -- folder where all processed data will be stored.  Example: '/data13/mkarlsso/Con/'
% FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.
% DAYS -- this is a vector containing the day numbers to add the information to.
% EPOCHS -- a vector containing the epoch numbers to add info to, example: [2 4 6]
% Properties:
%   These items are added as fields to the structure.  Two reserved
%   properties are:
%   'exposure' - must be the same length as the number of days * number of epochs (counts number of epoch exposures on the environment)
%   'exposureday' - must be the same length as the number of days
%   'experimentday' - must be the same length as the number of days
%   Any other names that the user wants to add are processed and added to the
%   structure.  The value will be copied across all epochs.
%
% Example: addtaskinfo('/data19/mkarlsso/Ale/','ale',1:10,[2 4],'exposures',1:20,'exposuredays',1:10,'description','TrackB')

exposure = [];
exposureday = [];
fieldvar = [];
fieldvarcount = 0;

daycount = 0;
epochcount = 0;

for option = 1:2:length(varargin)-1
    fieldvarcount = fieldvarcount + 1;
    if isstr(varargin{option})
        fieldvar(fieldvarcount).name = varargin{option};
        fieldvar(fieldvarcount).value = varargin{option+1};
    else
        error('Fields must be strings');
    end
end

for i = 1:length(days)
    if (days(i) < 10)
        dstring = '0';
    else
        dstring = '';
    end
    
    filename = fullfile(animdirect,[fileprefix,'task',dstring,num2str(days(i))]);
    try
        load(filename);
    catch
        task = [];
    end
    for j = 1:length(epochs)
        epochcount = epochcount+1;
        if (length(task{days(i)}) < epochs(j))
            task{days(i)}{epochs(j)} = [];
        end
        if ~isempty(fieldvar)
            for k = 1:length(fieldvar)
                    switch fieldvar(k).name
                        case {'experimentday', 'exposureday'} %disperse along the days
                            if length(fieldvar(k).value) ~= (length(days))
                                error(['The length of the', fieldvar(k).name, 'vector must equal the number of days']);
                            end
                            task{days(i)}{epochs(j)} = setfield(task{days(i)}{epochs(j)},fieldvar(k).name,fieldvar(k).value(i));
                        case 'exposure' %disperse along the epochs
                            if length(fieldvar(k).value) ~= (length(days)*length(epochs))
                                error(['The length of the', fieldvar(k).name, 'vector must equal the number of days * the number of epochs']);
                            end
                            task{days(i)}{epochs(j)} = setfield(task{days(i)}{epochs(j)},fieldvar(k).name,fieldvar(k).value(epochcount));
                        case 'dailyexposure'
                            if length(fieldvar(k).value) ~= (length(epochs))
                                error(['The length of the', fieldvar(k).name, 'vector must equal the number of epochs']);
                            end
                            task{days(i)}{epochs(j)} = setfield(task{days(i)}{epochs(j)},fieldvar(k).name,fieldvar(k).value(j));
                        otherwise %copy to every epoch
                            task{days(i)}{epochs(j)} = setfield(task{days(i)}{epochs(j)},fieldvar(k).name,fieldvar(k).value);    
                    end
            end
        end
    end
            
    save(filename,'task');
    
end
      

