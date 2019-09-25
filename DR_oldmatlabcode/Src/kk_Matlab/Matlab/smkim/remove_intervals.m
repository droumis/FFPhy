function new_intervals = remove_intervals(intervals,func)
%REMOVE_INTERVALS Given a set of real cadlag intervals, remove those that meet a test condition
%
%   NEW_INTERVALS = REMOVE_INTERVALS(INTERVALS,FUNC) removes the rows of
%   INTERVALS that satisfy a test condition FUNC. FUNC must be a function
%   handle which accepts an Nx2 real numeric array and returns an N-element
%   logical vector.
%
%Depends on:
%   IS_INTERVALS (written by smk)
%
%Written by smk, 2009 November 20.
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

if isempty(intervals)
  new_intervals = intervals;
  return;
end

if ~isa(func,'function_handle')
  error('FUNC must be a function handle');
end
try
  tf = func(intervals);
  assert(islogical(tf) && isvector(tf) && numel(tf) == size(intervals,1));
catch
  error('FUNC does not evaluate to a logical vector of the expected size');
end

new_intervals = intervals(find(~tf),:);

