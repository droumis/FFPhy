function new_intervals = join_intervals(intervals,threshold)
%JOIN_INTERVALS Consolidate intervals that are spaced apart less than a threshold
%
%   JOIN_INTERVALS(INTERVALS,THRESHOLD) takes the set of real cadlag intervals
%   represented by INTERVALS and consolidates neighboring intervals (rows) that
%   are less than THRESHOLD apart. THRESHOLD must be a real finite positive
%   scalar of the same numeric class as INTERVALS.
%
%Depends on:
%   IS_INTERVALS
%
%Written by SMK, 2009 November 20.
%

if (exist('is_intervals') ~= 2)
  error(['JOIN_INTERVALS depends on m-file IS_INTERVALS ' ...
      '(written by smk)']);
end

if ~is_intervals(intervals) 
  error(['INTERVALS argument does not appear to be a valid set of ' ...
      'non-overlapping real cadlag intervals']);
end
numeric_class = class(intervals);

if ~isnumeric(threshold) || ~isreal(threshold) || ~isscalar(threshold) || ...
    ~(threshold > 0) || ~isfinite(threshold)
  error('THRESHOLD argument must be real numeric positive finite scalar');
end
if ~isa(threshold,numeric_class)
  error('THRESHOLD and INTERVALS must be the same numeric class');
end

if (size(intervals,1) < 2)
  new_intervals = intervals;
  return;
end

gaps = intervals(2:end,1) - intervals(1:end-1,2);
% sanity check
assert(all(gaps > 0));
% Which gaps between rows of threshold are to be joined?
join_idx = (gaps < threshold);
new_intervals = zeros([size(intervals,1)-nnz(join_idx), 2],numeric_class);
% How do the rows of intervals map to the rows of new_intervals?
map_idx = (1:size(intervals,1))' - [0; cumsum(join_idx)];

for i = 1:size(new_intervals,1)
  new_intervals(i,1) = min(intervals(find(map_idx == i),1));
  new_intervals(i,2) = max(intervals(find(map_idx == i),2));
end





