function resettaskinfo(animdirect, fileprefix, days, keepers)
% function resettaskinfo(animdirect, fileprefix, days)
% RESETTASKINFO(animdirect, fileprefix, days, runepochs)
% function resettaskinfo(animdirect, fileprefix, days, keepers)

% If 'keepers' is not included, reset task structure to contain
% only 'type' and 'linearcoord' fields.
% 
% Otherwise, keep fields in cellarray 'keepers'
%
% ANIMDIRECT -- folder where all processed data will be stored.  Example: '/data13/mkarlsso/Con/'
% FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.
% DAYS -- this is a vector containing the day numbers to add the information to.

if nargin == 3
  keepers = {'type', 'linearcoord'};
elseif nargin ~= 4
  error('resettaskinfo: Wrong number of arguments');
end

for i = 1:length(days)
    
    filename = fullfile(animdirect,[fileprefix,'task',sprintf('%02d.mat',days(i))]);  
    if exist(filename,'file')
      load(filename);
    else
      fprintf('file not found: %s\n',filename);
      continue;
    end
    
    for j = 1:length(task{i})
        if isstruct(task{i}{j})
            names = fieldnames(task{i}{j});
            for k = 1:length(names)
                switch names{k}
                    case keepers

                    otherwise
                        task{i}{j} = rmfield(task{i}{j},names{k});
                end
            end
        end
    end
            
    save(filename,'task');  
end
      
