function excludePeriods = getExcludePeriods(time, included)
% excludePeriods = getExcludePeriods(time, included)
% Calculates the start and end times for all exclusion periods, given a
% time vector and an include vector of 1s and 0s of the same length.

if (length(time) ~= length(included))
    error('The TIME and INCLUDED vectors must me the same length');
end

starttimes = find((diff(included) == -1))+1;
starttimes = starttimes(:);
endtimes = find((diff(included) == 1));
endtimes = endtimes(:);
if (included(1) == 0)
    starttimes = [1; starttimes];
end
if (included(end) == 0)
    endtimes = [endtimes; length(included)];
end

excludePeriods = [time(starttimes) time(endtimes)];