function intervals = intersect_intervals(varargin)
%INTERSECT_INTERVALS Intersection of real cadlag interval sets.
%
%   INTERSECT_INTERVALS(A,B), for real cadlag interval sets A and B, returns
%   the set of intervals that are included in both A and B. An error is raised
%   if A and B are of different numeric classes.
%
%   INTERSECT_INTERVALS(A,B,C,...) is the same as 
%   INTERSECT_INTERVALS(A,INTERSECT_INTERVALS(B,INTERSECT_INTERVALS(C,...)))
%
%   See also IS_INTERVALS, FIND_INTERVALS, ISMEMBER_INTERVALS, DIFF_INTERVALS,
%   UNION_INTERVALS, XOR_INTERVALS written by smk.
%
%Depends on:
%   IS_INTERVALS (written by smk)
%
%Written by smk, 2009 June 23.
%

if (exist('is_intervals') ~= 2)
  error(['INTERSECT_INTERVALS depends on m-file IS_INTERVALS ' ...
      '(written by smk)']);
end

if (length(varargin) < 2)
  error('INTERSECT_INTERVALS takes at least two arguments');
end
if ~all(cellfun(@is_intervals,varargin))
  error(['At least one of the input arguments does not appear to be a ' ...
      'valid set of non-overlapping timestamp intervals']);
end
numeric_class = class(varargin{1});
if ~all(cellfun(@(arg) isa(arg,numeric_class),varargin))
  error('Arguments must be of the same numeric class');
end

% Recursive function composition
while (length(varargin) > 1)
  a = varargin{end-1};
  b = varargin{end};
  if isempty(a) || isempty(b)
    intervals = zeros([0 2],numeric_class);
  else
    % concatenate and sort
    x = sortrows([a; b]);
    assert(all(diff(x(:,1)) >= 0)); % make sure that we sorted correctly!

    % pre-allocate an output vector with adequate capacity
    intervals = zeros([size(a,1)+size(b,1) 2],numeric_class);
    nrows = 0; % running count of how many entries 

    i = 1;
    while (i < size(x,1))
      % Given interval x(i,:), we can find all overlapping intervals x(j,:) by
      % exploiting the property that x(:,1) is sorted
      j = i + find(x((i+1):end,1) < x(i,2));
      % Note that any such overlaps must occur between an interval in A and an
      % interval in B; intervals that belong to A do not overlap with each other,
      % and intervals that belong to B do not overlap with each other. Because
      % there are only two sets, A and B, if interval_1 overlaps with interval_2,
      % and interval_2 overlaps with interval_3, then we are guaranteed that
      % interval_1 does not overlap with interval_3.
      intervals((nrows+1):(nrows+numel(j)),:) = ...
          [bsxfun(@max,x(j,1),x(i,1)) bsxfun(@min,x(j,2),x(i,2))];
      nrows = nrows + numel(j);
      i = i + 1;
    end
    % trim excess
    intervals = intervals(1:nrows,:);
    if ~is_intervals(intervals)
      error('There is a bug in either INTERSECT_INTERVALS or IS_INTERVALS');
    end
  end
  % Pop last cell of varargin and replace the next cell in the stack with
  % the current value
  varargin(end) = [];
  varargin{end} = intervals;
end

