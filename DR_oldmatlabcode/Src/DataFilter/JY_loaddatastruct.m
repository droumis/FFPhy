function out = JY_loaddatastruct(animaldir, animalprefix, datatype, days)
% out = loaddatastruct(animaldir, animalprefix, datatype)
% out = loaddatastruct(animaldir, animalprefix, datatype, days)
%
% Load the components of a data cell array of combines then into one
% variable.  Datatype is a string with the name of the data structure in the loaded file.
% If DAYS is omitted, all files are loaded.  Otherwise, only the files for the
% specified days will be included.
% 20111221: changed file opening lines

if (nargin < 4)
    days = [];
end
out = [];
datafiles = dir([animaldir, animalprefix,datatype,'*.mat']); % removed 'datatype' since task file in reroute task is not stored there

for i = 1:length(datafiles)
    if isempty(days)
        s = datafiles(i).name;
        if isempty(strfind(s,'task')) && isempty(strfind(s,'linpos')); % do not load linpos or task files
            
            load([animaldir,datafiles(i).name]);
            eval(['out = datavaradd(out,',datatype,');']);
        end
        
    else
        s = datafiles(i).name;
        if isempty(strfind(s,'task')) && isempty(strfind(s,'linpos')); % do not load linpos or task files
            fileday = str2num(s(strfind(s,'_')+1:strfind(s,'.')-1));  %get the experiment day from the filename
            if (isempty(fileday))|(ismember(fileday,days))
                
                load([animaldir,datafiles(i).name]);
                eval(['out = datavaradd(out,',datatype,');']);
            end
        end
    end
end


%--------------------------------------
function out = datavaradd(origvar, addcell)

out = origvar;
for i = 1:length(addcell)
    if (~isempty(addcell{i}))
        out{i} = addcell{i};
    end
end

