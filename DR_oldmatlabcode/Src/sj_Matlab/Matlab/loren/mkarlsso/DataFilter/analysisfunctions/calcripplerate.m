function out = calcripplerate(index, excludetimes, pos, ripples, cellinfo, varargin)
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

if ~isempty(excludetimes)
    totaltime = (pos{index(1)}{index(2)}.data(end,1)-pos{index(1)}{index(2)}.data(1,1)) - sum(excludetimes(:,2)-excludetimes(:,1));
else
    totaltime = (pos{index(1)}{index(2)}.data(end,1)-pos{index(1)}{index(2)}.data(1,1));
end
riptimes = getripples(index, ripples, cellinfo, 'cellfilter', '(isequal($area, ''CA1'') | isequal($area, ''CA3''))','excludeperiods', excludetimes,'minstd',3);


rate = size(riptimes,1)/totaltime;

if (appendindex)
    out = [index rate]; %append the cell index to the value
else
    out = rate;
end