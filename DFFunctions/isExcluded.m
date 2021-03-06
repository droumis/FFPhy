function out = isExcluded(times, excludePeriods)
% out = isExcluded(times, excludePeriods)
% Returns a list of 1s and 0s signifying inclusion (0) or exclusion (1), as defined
% by the exclude start and end times in excludePeriods.  ExcludePeriods is
% an n by 2 matrix where the columns are the start and end times for each
% exclusion period.
% 1 signifies times that fall in between start and end times of excludePeriods
% so if you want indices into times that are not excluded, you want to do
% ~isExcluded
if ~isempty(excludePeriods)
    oneborder = [(excludePeriods(:,1)-.0000001);excludePeriods(:,2)+.0000001];
    oneborder(:,2) = 0;
    zeroborder = excludePeriods(:);
    zeroborder(:,2) = 1;
    sortedMatrix = [[-inf 0]; sortrows([oneborder;zeroborder],1); [inf 0]];
else
    sortedMatrix = [[-inf 0]; [inf 0]];
end
out = sortedMatrix(lookup(times,sortedMatrix(:,1)),2);

