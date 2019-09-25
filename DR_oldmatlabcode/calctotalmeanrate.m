function out = calctotalmeanrate(indices, excludetimes, spikes, varargin)
% out = calctotalmeanrate(index, excludetimes, spikes, options)
% Calculates the mean rate of a cell for a given epoch as the number of
% spikes/ time.  Excluded time periods are not included in the total time.
%
% Options:
%   'appendindex', 1 or 0 -- set to 1 to append the cell index to the
%   output [tetrode cell value].  Default 0.
%



out = [];  
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


for cellNum = 1:size(indices.spikes,1)
    
    index = [indices.epochs indices.spikes(cellNum,:)];
    
    
    timeRange = spikes{index(1)}{index(2)}{index(3)}{index(4)}.timerange./10000; %units in 10000 seconds
    ontime = diff(timeRange);
    if ~isempty(excludetimes)
        
        excludetimes(find(excludetimes > timeRange(2))) = timeRange(2);
        excludetimes(find(excludetimes < timeRange(1))) = timeRange(1);
        totalexclude = sum(excludetimes(:,2) - excludetimes(:,1));
    else
        totalexclude = 0;
    end
    
    totalontime = ontime-totalexclude;
    
    
    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
        if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.time)
            %Use 'isExcluded' to determine which spikes are excluded in the excludetimes list 
            goodspikes = ~isExcluded(spikes{index(1)}{index(2)}{index(3)}{index(4)}.time(:,1), excludetimes);
        else
            goodspikes = 0;
        end
    else
        goodspikes = nan;
    end
    numgoodspikes = sum(goodspikes);
    
    if (totalontime > 1)
        rate = numgoodspikes/totalontime;
    else
        rate = 0;
    end
    if (appendindex)
        out = [out;[index rate]]; %append the cell index to the value
    else
        out = [out; rate];
    end
end


