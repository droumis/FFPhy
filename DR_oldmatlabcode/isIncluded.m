function out = isIncluded(times, includePeriods)
% Similar to isExcluded - SJ. Dec 2013. For use in DFAsj_getthetacrosscov_timecondition
% out = isExcluded(times, excludePeriods)

% Returns a list of 1s and 0s signifying inclusion or exclusion, as defined
% by the include start and end times in includePeriods.  IncludePeriods is
% an n by 2 matrix where the columns are the start and end times for each
% inclusion period.
% 1 signifies times that fall in between start and end times of
% includePeriods

if ~isempty(includePeriods)
    oneborder = [(includePeriods(:,1)-.0000001);includePeriods(:,2)+.0000001];
    oneborder(:,2) = 0;
    zeroborder = includePeriods(:);
    zeroborder(:,2) = 1;
    sortedMatrix = [[-inf 0]; sortrows([oneborder;zeroborder],1); [inf 0]];
else
    sortedMatrix = [[-inf 0]; [inf 0]];
end
out = sortedMatrix(lookup(times,sortedMatrix(:,1)),2);

