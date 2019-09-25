function [bool, loc] = ismember_intervals(x,intervals)
%ISMEMBER_INTERVALS Test whether array elements fall within real cadlag interval set.
%
%   BOOL = ISMEMBER_INTERVALS(X,INTERVALS), for real array X, returns true for
%   elements of X that lie within the real cadlag intervals represented by
%   INTERVALS. An error is raised if INTERVALS and X are not of the same
%   numeric class. Note that the time intervals are assumed to be left-closed,
%   right-open (cadlag) intervals; for the ith row of INTERVALS, the test
%   condition is
%       (X >= INTERVALS(i,1)) & (X < INTERVALS(i,2))
%   
%   The return value BOOL is a logical array of the same size as X.
%
%   [BOOL, LOC] = ISMEMBER_INTERVALS(X,INTERVALS) also returns an array LOC
%   which contains, for each element of X, the subscript of the row of
%   INTERVALS which bounds the element or 0 if there is no such row.
%
%   See also IS_INTERVALS, FIND_INTERVALS, DIFF_INTERVALS, INTERSECT_INTERVALS,
%   UNION_INTERVALS, XOR_INTERVALS written by smk.
%
%Depends on:
%   IS_INTERVALS (written by smk)
%
%Written by smk, 2009 June 23.
%

if (exist('is_intervals') ~= 2)
  error(['ISMEMBER_INTERVALS depends on m-file IS_INTERVALS ' ...
      '(written by smk)']);
end

if ~isnumeric(x) || ~isreal(x)
  error('X argument must be real numeric array');
end
if ~is_intervals(intervals) 
  error(['INTERVALS argument does not appear to be a valid set of ' ...
      'non-overlapping real cadlag intervals']);
end
numeric_class = class(intervals);
if ~isa(x,numeric_class)
  error('X and INTERVALS must be the same numeric class');
end

% Allocate output of the same size as X
bool = false(size(x));

if isempty(intervals) || isempty(x)
  return;
end

% Combine all elements of X and INTERVALS into a single column vector.
% MATLAB documentation for SORT function says:
%
%    "When more than one element has the same value, the order of the
%    elements are preserved in the sorted result and the indexes of
%    equal elements will be ascending in any index matrix."
%
% The order of concatenation is very important!
t = [ intervals(:,1); intervals(:,2); x(:) ];
% Matching vector of indicator values. +Inf indicates the starting edge of an
% interval, and -Inf indicates the ending edge of an interval.
i = [ +Inf([size(intervals,1) 1]); -Inf([size(intervals,1) 1]); ...
    (1:numel(x))' ];
[t, sortorder] = sort(t);
i = i(sortorder);

% This last part is subtle. Break into steps and study a toy example if you
% don't understand
mask = isfinite(i);
bool(i(cumsum((i == +Inf) - (i == -Inf)) & mask)) = true;

if nargout == 2
  loc = zeros(size(x));
  % Which time interval was it in? Count using CUMSUM.
  n = ceil(cumsum(~mask)/2);
  loc(i(mask)) = n(mask);
  % But wait! The CUMSUM count erroneously assigns non-zero indices to elements
  % of X which do not fall within any time interval. Use the first
  % output argument to replace these elements with zeros.
  loc(~bool) = 0;
  assert(isequal(loc > 0,bool));
  assert(all(ismember(unique(loc),0:size(intervals,1))));
end


