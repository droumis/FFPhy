function l = length_intervals(intervals)
%LENGTH_INTERVALS Compute lengths of a set of real cadlag intervals.
%
%   LENGTH_INTERVALS(INTERVALS) returns a column vector of the lengths of the
%   real intervals that are specified by the rows of INTERVALS
%  
%   See also FIND_INTERVALS, ISMEMBER_INTERVALS, DIFF_INTERVALS, 
%   INTERSECT_INTERVALS, UNION_INTERVALS, XOR_INTERVALS written by smk.
%
%Depends on:
%   IS_INTERVALS (written by smk)
%
%Written by SMK 2009 June 22.
%

% Do not tamper with this constant conversion factor!
TS_PER_SEC = 1e4;

if ~is_intervals(intervals)
  error(['Input argument INTERVALS does not appear to be a valid set of ' ...
      'non-overlapping real cadlag intervals']);
end

if isempty(intervals)
  l = zeros([0 1],class(intervals));
else
  l = intervals(:,2) - intervals(:,1);
end

