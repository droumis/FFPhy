function out = isIncluded(times, intervals, varargin)
% out = isIncluded(times, intervals)
% @DKR.. correction of isExcluded. bc wtf was someone thinking when they
% named that.. 
% Returns a list of 1s and 0s signifying inclusion or exclusion, as defined
% by the interval start and end times in intervals.  Intervals is
% an n by 2 matrix where the columns are the start and end times for each
% interval.
% 1 signifies times that fall in between start and end times of intervals

offset = .002; % seconds buffer
if ~isempty(varargin)
    assign(varargin{:})
end

if ~isempty(intervals)
    oneborder = [(intervals(:,1)-offset); intervals(:,2)+offset];
    oneborder(:,2) = 0;
    zeroborder = intervals(:);
    zeroborder(:,2) = 1;
    sortedMatrix = [[-inf 0]; sortrows([oneborder;zeroborder],1); [inf 0]];
else
    sortedMatrix = [[-inf 0]; [inf 0]];
end
out = logical(sortedMatrix(lookup(times,sortedMatrix(:,1)),2));

