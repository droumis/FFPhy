function [out] = nInInterval(times, excludePeriods, incends)
%out = nInInterval(times, excludePeriods, includeStartEnd)
% Returns a list of numbers, one per time interval bracketed by exlude periods,
% of the number of times in each included period.
%
% ExcludePeriods is an n by 2 matrix where the columns are the start and end
% times for each exclusion period.
%
% includeStartEnd is 0 if the exclude periods begin and end at the beginning
% and end of the recording period, in which case we don't include the times
% before the first exclude time or after the last time. A value of 1 indicates
% that we should include those times.

if ~isempty(excludePeriods)
    incborder = sort([(excludePeriods(:,1)-.0000001);excludePeriods(:,2)+.0000001]);
    incborder(:,2) = 3:(size(incborder,1)+2);
    zeroborder = excludePeriods(:);
    zeroborder(:,2) = 0;
    sortedMatrix = [[-inf 2]; sortrows([incborder;zeroborder],1); [inf ...
        (incborder(end,2)+1)]];
    if (incends == 0)
        % get rid of the first two and last two elements
        sortedMatrix = sortedMatrix(3:(end-2),:);
    end
    out = sortedMatrix(lookup(times,sortedMatrix(:,1)),2);
    % the result of the convoluted procedure is that each time in an included
    % interval will be between 2n and 2n+1, so we can figure out which interval
    % it is by taking the floor of the value / 2
    interval = floor(out/2);
    % now we histogram interval from 1 to the number of included periods
    % there is one more included than excluded period
    nintervals = size(excludePeriods,1) + 1;
    %if ~isempty(interval) & sum(interval)~=0
        [out bins] = flhist(interval, 0.5:1:(nintervals+0.5));
    %end
else
    out = length(times);
end
