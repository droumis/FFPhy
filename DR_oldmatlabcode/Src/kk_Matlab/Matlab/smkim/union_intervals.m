function intervals = union_intervals(varargin)
%UNION_INTERVALS Union of real cadlag interval sets.
%
%   UNION_INTERVALS(A,B), for real cadlag interval sets A and B, returns the
%   set of intervals that cover all intervals in either A or B or both. An
%   error is raised if A and B are of different numeric classes.
%
%   UNION_INTERVALS(A,B,C,...) is the same as 
%   UNION_INTERVALS(A,UNION_INTERVALS(B,UNION_INTERVALS(C,...)))
%
%   See also IS_INTERVALS, FIND_INTERVALS, ISMEMBER_INTERVALS, DIFF_INTERVALS,
%   INTERSECT_INTERVALS, XOR_INTERVALS written by smk.
%
%Depends on:
%   IS_INTERVALS (written by smk)
%
%Written by smk, 2009 June 23.
%

if (exist('is_intervals') ~= 2)
  error(['UNION_INTERVALS depends on m-file IS_INTERVALS ' ...
      '(written by smk)']);
end

if (length(varargin) < 1)
  error('UNION_INTERVALS takes at least two arguments');
end
if ~all(cellfun(@is_intervals,varargin))
  error(['At least one of the input arguments does not appear to be a ' ...
      'valid set of non-overlapping timestamp intervals']);
end
numeric_class = class(varargin{1});
if ~all(cellfun(@(arg) isa(arg,numeric_class),varargin))
  error('Arguments must be of the same numeric class');
end

if (length(varargin) == 1)
  intervals = varargin{1};
end

% Recursive function composition
while (length(varargin) > 1)
  a = varargin{end-1};
  b = varargin{end};
  if isempty(a)
    intervals = b;
  elseif isempty(b)
    intervals = a;
  else
    % concatenate and sort
    intervals = sortrows([a; b]);
    assert(all(diff(intervals(:,1)) >= 0)); % make sure that we sorted correctly!
    i = 1;
    while (i < size(intervals,1))
      % We can exploit the property that intervals(:,1) is sorted to perform
      % iterative greedy consolidation
      if intervals(i,2) >= intervals(i+1,1)
        intervals(i,2) = max(intervals(i:(i+1),2));
        intervals(i+1,:) = []; % don't increment i after deleting a row!
      else
        i = i + 1;
      end
    end
    if ~is_intervals(intervals)
      error('There is a bug in either UNION_INTERVALS or IS_INTERVALS');
    end
  end
  % Pop last cell of varargin and replace the next cell in the stack with
  % the current value
  varargin(end) = [];
  varargin{end} = intervals;
end

