function bool = is_intervals(intervals)
%IS_INTERVALS Validates whether the input is a valid set of non-overlapping real cadlag intervals
%
%   IS_INTERVALS(INTERVALS) returns true if INTERVALS is an Nx2 (N >= 0) real
%   numeric array, where each row corresponds to an interval 
%   [start_value, end_value] such that (end_value > start_value), all intervals
%   are non-overlapping, the rows are in monotonically-increasing order, and no
%   elements are NaN or infinite. The intervals are assumed to be left-closed,
%   right-open (cadlag) intervals. The return value is a logical scalar.
%  
%   See also FIND_INTERVALS, ISMEMBER_INTERVALS, DIFF_INTERVALS, 
%   INTERSECT_INTERVALS, UNION_INTERVALS, XOR_INTERVALS written by smk.
% 
%Written by SMK 2009 June 22.
%

if ~isnumeric(intervals) || ~isreal(intervals) || ...
    (size(intervals,2) ~= 2) || ...
    ~all(intervals(:,1) < intervals(:,2)) || ...
    ((size(intervals,1) > 1) && ...
    ~all(intervals(2:end,1) >= intervals(1:end-1,2))) || ...
    any(~isfinite(intervals(:))) || any(isnan(intervals(:)))
  bool = false;

else
  bool = true;

end

