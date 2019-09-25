function out = JY_loaddatastructlinpos(animaldir, animalprefix, datatype, daysind)
% out = loaddatastruct(animaldir, animalprefix, datatype)
% out = loaddatastruct(animaldir, animalprefix, datatype, days)
%
% Load the components of a data cell array in a linpos file, combines then into one
% variable.  Datatype is a string with the name of the data structure in the loaded file.  
% If DAYS is omitted, all files are loaded.  Otherwise, only the files for the
% specified days will be included.
% 20111221: changed file opening lines

if (nargin < 4)
    daysind = [];
end
out = [];
datafiles = dir([animaldir, animalprefix,'*linpos.mat']); % removed 'datatype' since task file in reroute task is not stored there

for i = 1:length(datafiles)
    if isempty(daysind)
        load([animaldir,datafiles(i).name]);
        eval(['out = datavaradd(out,',datatype,');']);
        
    else
        s = datafiles(i).name;
        
        fileday = str2num(s(strfind(s,'_')+1:strfind(s,'.')-7));  %get the experiment day from the filename
        if (isempty(fileday)) | (ismember(fileday,daysind))
            load([animaldir,datafiles(i).name]);
            eval(['out = datavaradd(out,',datatype,');']);
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
        
