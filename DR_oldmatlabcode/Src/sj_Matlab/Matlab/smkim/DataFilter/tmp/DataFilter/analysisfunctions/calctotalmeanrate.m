function out = calctotalmeanrate(index, excludetimes, spikes, varargin)
% out = calctotalmeanrate(index, excludetimes, spikes, options)
% Calculates the mean rate of a cell for a given epoch as the number of
% spikes/ time.  Excluded time periods are not included in the total time.
%
% Options:
%   'appendindex', 1 or 0 -- set to 1 to append the cell index to the
%   output [tetrode cell value].  Default 0.
%


appendindex = 0;
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

ontime = diff(spikes{index(1)}{index(2)}{index(3)}{index(4)}.timerange)/10000;
if ~isempty(excludetimes)
    totalexclude = sum(excludetimes(:,2) - excludetimes(:,1));
else
    totalexclude = 0;
end

totalontime = ontime-totalexclude;
if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
        goodspikes = ~isExcluded(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1), excludetimes);
    else
        goodspikes = 0;
    end
else
    goodspikes = nan;
end
numgoodspikes = sum(goodspikes);
rate = numgoodspikes/totalontime;

if (appendindex)
    out = [index rate]; %append the cell index to the value
else
    out = rate;
end