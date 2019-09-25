function out = loadeegstruct(animaldir, animalprefix, datatype, days, epochs, tet)

% out = loadeegstruct(animaldir, animalprefix, datatype, days, epochs, tetrodes)
%
% Load the components of an eeg cell array of combines then into one
% variable. 
%	datatype is a string with the base name of the files (e.g. 'theta') 
% 	days specifies the list of days to load 
% 	epochs specifies the list of epochs to load 
% 	tetrodes specifies the list of tetrodes to load 
%
%	if epochs is empty, all epochs (1-10) are loaded
%	if tetrodes is empty, all tetrodes (1-50) are loaded
%
% Be aware that eeg files tend to be large, so loading a large number of them
% is likely to lead to problems

if ~exist('epochs') || isempty(epochs) % || short circuits!!!
    epochs = 1:50;
end
if ~exist('tet') || isempty(tet) % || short circuits!!!  
    tet = 1:50;
end

if size(days,1) > 1
   days = days';
end

if size(epochs,1) > 1
   epochs = epochs';
end

if size(tet,1) > 1
   tet = tet';
end

out = [];
% create the list of data files 
for d = days
     for e = epochs
	for t = tet
       if ~ischar(animalprefix) | ~ischar(datatype)
          keyboard
       end
       if length(d) > 1 | length(e) > 1 | length(t) > 1
          keyboard
       end
	    fname = fullfile(animaldir,'EEG',sprintf('%s%s%02d-%d-%02d.mat', ...
	    		animalprefix, datatype, d, e, t));
       if exist(fname,'file');
          try
            load(fname);
            eval(['out = datavaradd(out,',datatype,');']);
          catch
            warning('load failed for some reason');
          end
       end
   end
     end
end


%--------------------------------------
function out = datavaradd(origvar, addcell)

out = origvar;
for d = 1:length(addcell)
    for e = 1:length(addcell{d})
	for t = 1:length(addcell{d}{e})
	    if (~isempty(addcell{d}{e}{t}))
		out{d}{e}{t} = addcell{d}{e}{t};
	    end
	end
    end
end
