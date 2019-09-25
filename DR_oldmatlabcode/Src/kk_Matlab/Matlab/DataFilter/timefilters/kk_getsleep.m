function [out] = kk_getsleep(animaldir,animalprefix, epochs, tetlist, varargin)
% out = getriptimes(animaldir,animalprefix,epochs, tetlist, options)

%%    modified by kk May 2013 to avoid end ripple time problems
%
%     animaldir and animal prefix are strings indicating the base director for
%     the animal's data and the prefix for the data files
%
%     epochs is an Nx2 list of days and epochs
%
%     tetlist is a list of tetrodes to use or an empty matrix if the
%     'cellfilter' option is used.
% Produces a cell structure with a time field and an nripples field which
% indicates the number of electrodes with a ripple at each time point
%
% Examples:
% getriptimes('/data/name/Fre', 'fre', epochs, 1)
% getriptimes('/data/name/Fre', 'fre', epochs, [], 'cellfilter', '(isequal($area, ''CA1''))')

% assign the options

for option = 1:2:length(varargin)-1
    switch varargin{option}
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

loaddays = unique(epochs(:,1));
sleep = loaddatastruct(animaldir, animalprefix, 'sleep', loaddays);

for i = 1:size(epochs,1)
    
    % setup vector of 1 ms points
	times = sleep.timerange(1):0.001:sleep.timerange(end);
    % retrieve sleep times
    sleepperiods = [ sleep{epochs(i,1)}{epochs(i,2)}.starttime   sleep{epochs(i,1)}{epochs(i,2)}.endtime ];
    % output
    out{epochs(i,1)}{epochs(i,2)}.time = times;
    out{epochs(i,1)}{epochs(i,2)}.sleep = list2vec(sleepperiods,times);
    
    clear times
  
end
