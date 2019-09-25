function intervals = xor_intervals(a,b)
%XOR_INTERVALS Exclusive-or of cadlag real interval sets.
%
%   XOR_INTERVALS(A,B), for cadlag real intervals sets A and B, returns the set
%   of intervals that cover all times that are in either A or B, but not in
%   both. An error is raised if A and B are of different numeric classes.
%
%   See also IS_INTERVALS, FIND_INTERVALS, ISMEMBER_INTERVALS, DIFF_INTERVALS,
%   INTERSECT_INTERVALS, UNION_INTERVALS written by smk.
%
%Depends on:
%   IS_INTERVALS (written by smk)
%   UNION_INTERVALS (written by smk)
%   INTERSECT_INTERVALS (written by smk)
%   DIFF_INTERVALS (written by smk)
%
%Written by smk, 2009 June 23.
%

if (exist('is_intervals') ~= 2)
  error('XOR_INTERVALS depends on m-file IS_INTERVALS (written by smk)');
end
if (exist('union_intervals') ~= 2)
  error('XOR_INTERVALS depends on m-file UNION_INTERVALS (written by smk)');
end
if (exist('intersect_intervals') ~= 2)
  error('XOR_INTERVALS depends on m-file INTERSECT_INTERVALS (written by smk)');
end
if (exist('diff_intervals') ~= 2)
  error('XOR_INTERVALS depends on m-file DIFF_INTERVALS (written by smk)');
end

if ~is_intervals(a) || ~is_intervals(b)
  error(['At least one of the input arguments does not appear to be a ' ...
      'valid set of non-overlapping timestamp intervals']);
end
numeric_class = class(a);
if ~isa(b,numeric_class)
  error('A and B must be of the same numeric class');
end

intervals = diff_intervals(union_intervals(a,b),intersect_intervals(a,b));

if ~is_intervals(intervals)
  error('There is a bug in either XOR_INTERVALS or IS_INTERVALS');
end



