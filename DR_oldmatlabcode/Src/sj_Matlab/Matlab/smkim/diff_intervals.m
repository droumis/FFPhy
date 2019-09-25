function intervals = diff_intervals(a,b)
%DIFF_INTERVALS Difference of real cadlag interval sets.
%
%   DIFF_INTERVALS(A,B), for real cadlag interval sets A and B, returns the set
%   of intervals that cover times that are in A but not in B. An error is
%   raised if A and B are of different numeric classes.
%
%   See also IS_INTERVALS, FIND_INTERVALS, ISMEMBER_INTERVALS, 
%   INTERSECT_INTERVALS, UNION_INTERVALS, XOR_INTERVALS written by smk.
%
%Depends on:
%   IS_INTERVALS (written by smk)
%   INTERSECT_INTERVALS (written by smk)
%
%Written by smk, 2009 June 23.
%

if (exist('is_intervals') ~= 2)
  error(['DIFF_INTERVALS depends on m-file IS_INTERVALS ' ...
      '(written by smk)']);
end
if (exist('intersect_intervals') ~= 2)
  error(['DIFF_INTERVALS depends on m-file INTERSECT_INTERVALS ' ...
      '(written by smk)']);
end

if ~is_intervals(a) || ~is_intervals(b)
  error(['At least one of the input arguments does not appear to be a ' ...
      'valid set of non-overlapping real intervals']);
end

numeric_class = class(a);
if ~isa(b,numeric_class)
  error('A and B must be of the same numeric class');
end

if isempty(b) || isempty(a)
  intervals = a;
  return;
end

% Construct the complement of B
c = zeros([size(b,1)+1 2],numeric_class);
if isfloat(c)
  c(1,1) = realmin(numeric_class);
  c(end,2) = realmax(numeric_class);
else
  c(1,1) = intmin(numeric_class);
  c(end,2) = intmax(numeric_class);
end
c(1:end-1,2) = b(:,1);
c(2:end,1) = b(:,2);
% Remove zero-length degenerate intervals
c(find(c(:,1) == c(:,2)),:) = [];
% Validate
if ~is_intervals(c)
  error('bug in code');
end

% Now take the intersection of A and (complement of B)
intervals = intersect_intervals(a,c);

if ~is_intervals(intervals)
  error('There is a bug in either DIFF_INTERVALS or IS_INTERVALS');
end

