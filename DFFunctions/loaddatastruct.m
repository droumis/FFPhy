function out = loaddatastruct(animaldir, animalprefix, datatype, days)

% out = loaddatastruct(animaldir, animalprefix, datatype)
% out = loaddatastruct(animaldir, animalprefix, datatype, days)
%
% Load the components of a data cell array of combines then into one
% variable.  Datatype is a string with the base name of the files.  If DAYS
% is omitted, all files are loaded.  Otherwise, only the files for the
% specified days will be included.

% instead of this function, use loadeegstruct to load eeg structures

if (nargin < 4)
    days = [];
end
out = [];
datafiles = dir([animaldir, animalprefix, datatype,'*']);
for i = 1:length(datafiles)
    if isempty(days)
        load([animaldir,datafiles(i).name]);
        %         try
        eval(['out = datavaradd(out,',datatype,');']);
        %         catch
        %             eval([datatype '= ' datatype 'info;']);
        %             eval(['out = datavaradd(out,',datatype,');']);
        %         end
    else
        s = datafiles(i).name;
        fileday = str2num(s(strfind(s,datatype)+length(datatype):strfind(s,'.')-1));  %get the experiment day from the filename
        if (isempty(fileday))|(ismember(fileday,days))
            load([animaldir,datafiles(i).name]);
            eval(['out = datavaradd(out,',datatype,');']);
        end
    end
end


%--------------------------------------
function out = datavaradd(origvar, addcell)

out = origvar;
try
    for i = 1:length(addcell)
        
        if (~isempty(addcell{i}))
            out{i} = addcell{i};
        end
    end
catch
    for i = 1:length(addcell.data)
        if (~isempty(addcell.data{i}))
            out{i} = addcell.data{i};
        end
    end
end

