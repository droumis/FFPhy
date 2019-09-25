function eegtet = geteegtet(adir, aprefix, ind, datatype, varargin)
%function eegtet = geteegtet(animaldir, animalprefix, ind, datatype, options)
% 
% Used as part of the eeg data filter, this function returns, for the given
% index and options, the eeg tetrode to use for analyzing the current tetrode
%
% animaldir is the directory for the animal's data
% animalprefix is the prefix name for the animal
% ind is the index of the cell being analyzed
% datatype is the datatype of the eeg file to be loaded (e.g. 'theta' or
% 	'gamma')
%
% options are
%	'sametet', 0 or 1	 
%		returns the tetrode number of the current cell 
%		default 1
%
%	'tetfilter', 'tetfilterstring'
%		applies tetfilter string to limit the set of tetrodes examined
%		
%	'maxvar', 0 or 1
%		returns the tetrode number of the tetrode with the maximal 
%		variance eeg
%		default 0
%	
%	'file', 'datafiletype' 
%		loads the file associated with datafiletype (e.g. 'thetaelect')
%               which should contain a variable called, for example, 
%               'eegtetrode' where eegtetrode{day}{epoch} specifies the 
%		tetrode to use.

% parse the options

maxvar = 0;
sametet = 1;
tetfilter = '';
file = '';

if (isempty(varargin))
    error('geteegtet requires at least on option to be specified');
end

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'sametet'
	    eegtet = ind(3);
	    return;
        case 'tetfilter'
	    tetfilter = varargin{option+1};
        case 'maxvar'
	    maxvar = varargin{option+1};
        case 'file'
	    try
		eegtetrode = loaddatastruct(adir, aprefix, varargin{option+1});
		eegtet = eegtetrode{ind(1)}{ind(2)};
	    catch
		error(sprintf('problem accessing eegtetrode{%d}{%d} from file %s', ind(1), ind(2), varargin{option+1}));
	    end
	    return;
        otherwise
            error(['Input ''', varargin{option}, ''' not defined']);
    end
end

% we will only get here if we need to do some calculations, so we get the set
% of tetrodes if tetfilter was specified and load up the list
if (~isempty(tetfilter))
    % this will cause us to ignore tetlist
    tetinfo = loaddatastruct(adir, aprefix, 'tetinfo');
    % Shantanu - CHange. This is wrong. Have to pass day and epoch
    tetlist = evaluatefilter(tetinfo{ind(1)}{ind(2)}, tetfilter);
else
    % set the tetlist to empty to load all tetrodes
    tetlist = [];
end


% calculate the specified measure to determine which channel to use
eegvar = zeros(length(tetlist),1);
if (maxvar)
    for i = 1:length(tetlist)
	eeg = loadeegstruct(adir, aprefix, datatype, ind(1), ind(2),tetlist(i,2));
	eegvar(i) = var(double(eeg{ind(1)}{ind(2)}{tetlist(i,2)}.data(:,1)));
    end
    [m, eind] = max(eegvar);
    eegtet = tetlist(eind,2);
    return;
else
    eegtet = tetlist;
end

