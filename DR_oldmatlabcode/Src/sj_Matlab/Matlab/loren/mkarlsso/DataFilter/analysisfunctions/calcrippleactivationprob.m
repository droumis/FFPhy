function activatedripples = calcrippleactivationprob(index, excludeperiods, spikes, ripples, cellinfo,varargin)
% activatedripples = calcrippleactivationprob(index, excludeperiods, spikes, ripples, cellinfo,varargin)
% Finds which ripples the cell fired at least one spike in.  Returns a
% vector of ones and zeros with the length of the number of ripples that
% occurred.
%
% Options:
% 'appendindex', 1: appends the cell index to the output
% 'ratio', 1: outputs the proportion of the ripples where the cell was
%  activated instead of the entire vector.
%

%Assign options
appendindex = 0;
ratioval = 0;
for option = 1:2:length(varargin)-1   
    if isstr(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'ratio'
                ratioval = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

%Find the ripples, using the tetrode in CA1 with the most cells, and a 3
%STD threshold
riptimes = getripples(index(1:2), ripples, cellinfo, 'cellfilter', 'isequal($area, ''CA1'')','excludeperiods', excludeperiods, 'minstd',3);
activatedripples = zeros(1,size(riptimes,1));

%Find the ripples in which the cell was activated
if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)})
    if ~isempty(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data)
        goodspikes = spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(find(~isExcluded(spikes{index(1)}{index(2)}{index(3)}{index(4)}.data(:,1), excludeperiods)),1);
        assignedripples = unique(periodAssign(goodspikes, riptimes(:,[1 3]))); 
        for i = 1:length(assignedripples)
            if (assignedripples(i) > 0)
                activatedripples(assignedripples(i)) = 1;
            end
        end
    end   
else
    activatedripples = [];
end


if ((ratioval) && ~isempty(activatedripples))
    %turn results into a ratio
    tmpl = length(activatedripples);
    tmpn = length(find(activatedripples));
    %activatedripples = [tmpn tmpl];
    activatedripples = tmpn/tmpl;
end
if (appendindex)
    activatedripples = [index activatedripples];
end



function out = periodAssign(times, periods)
% out = periodAssign(times, periods)
% TIMES is a vector of times
% PERIODS is an N by 2 list of start and end times
% Returns the index of the period that each time falls into.  If a time
% does not fall into one of the periods, a zero is returned.
% This function assumes that the periods are not overlapping.
%

if ~isempty(periods)
    oneborder = [(periods(:,1)-.0000001);periods(:,2)+.0000001];
    oneborder(:,2) = 0;
    insideborder = [(periods(:,1)+.0000001) (1:length(periods))'; (periods(:,2)-.0000001) (1:length(periods))'];    
    sortedMatrix = [[-inf 0]; sortrows([oneborder;insideborder],1); [inf 0]];
else
    sortedMatrix = [[-inf 0]; [inf 0]];
end
out = sortedMatrix(lookup(times,sortedMatrix(:,1)),2);
