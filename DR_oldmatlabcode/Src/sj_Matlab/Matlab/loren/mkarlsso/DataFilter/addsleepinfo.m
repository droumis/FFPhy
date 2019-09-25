function addsleepinfo(animdirect, fileprefix)

% ADDSLEEPINFO(animdirect, fileprefix)
%
% This function adds descriptive information to the 'task' structure of a
% given animal, specifically the sleep epochs. It saves the information for you.
%
% ANIMDIRECT -- folder where all processed data will be stored.  Example: '/data13/mkarlsso/Con/'
% FILEPREFIX -- also the first three letters of the animal's name (example 'con'), which will attach to
%              the beginning of the .mat files containing the variables.  I recommend using this.
%              If not, type ''.


currdir = pwd;
cd(animdirect);
files = dir([fileprefix,'task*']);
for i = 1:length(files)
    clear task
    load(files(i).name);
    for j = 1:length(task)
        sleepcount = 0;
        if (~isempty(task{j}))
            for k = 1:length(task{j})
                if (isfield(task{j}{k},'type') && (isequal(task{j}{k}.type,'sleep')))
                    sleepcount = sleepcount + 1;
                    task{j}{k}.sleepnum = sleepcount;
                    task{j}{k}.runbefore = [];
                    task{j}{k}.runafter = [];
                    if ((k > 1) && (isfield(task{j}{k-1},'type')) && (isequal(task{j}{k-1}.type,'run')))
                        task{j}{k}.runbefore = task{j}{k-1};
                    end
                    if ((k < length(task{j})) && (isfield(task{j}{k+1},'type')) && (isequal(task{j}{k+1}.type,'run')))
                        task{j}{k}.runafter = task{j}{k+1};
                    end
                    
                end
                
            end
        end
    end
    save(files(i).name,'task');
end




cd(currdir);