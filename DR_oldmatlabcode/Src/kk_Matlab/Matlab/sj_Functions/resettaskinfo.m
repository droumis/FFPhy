function resettaskinfo(animdirect, fileprefix, days)

% RESETTASKINFO(animdirect, fileprefix, days, runepochs)
%
% This function resets the task structure to only contain the 'type' and
% 'linearcoord' fields
%
% ANIMDIRECT -- folder where all processed data will be stored.  Example: '/data13/mkarlsso/Con/'
% FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.
% DAYS -- this is a vector containing the day numbers to add the information to.




for i = 1:length(days)
    if (days(i) < 10)
        dstring = '0';
    else
        dstring = '';
    end
    
    filename = [animdirect,fileprefix,'task',dstring,num2str(days(i))];  
    load(filename);
    
    for j = 1:length(task{i})
        if isstruct(task{i}{j})
            names = fieldnames(task{i}{j});
            for k = 1:length(names)
                switch names{k}
                    case {'type', 'linearcoord'}

                    otherwise
                        task{i}{j} = rmfield(task{i}{j},names{k});
                end
            end
        end
    end
            
    save(filename,'task');  
end
      