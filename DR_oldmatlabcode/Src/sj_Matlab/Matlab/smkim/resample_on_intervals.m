function [x_resampled, y_resampled] = resample_on_intervals(x,y,x_intervals,x_tol,varargin)
%RESAMPLE_ON_INTERVALS Resample function on non-overlapping intervals of domain.
%   [X_RESAMPLED, Y_RESAMPLED] = RESAMPLE_ON_INTERVALS(X,Y,X_INTERVALS,X_TOL)
%   takes a function Y which is sampled at points X in its domain and resamples
%   in on snippets that cover intervals X_INTERVALS. The resampling is done by
%   interpolation on a grid with grid spacing less than or equal to X_TOL. X
%   must be X must be a column vector of strictly monotonically increasing
%   unique real floating-point values, such that X(1) == min(X) and X(end) ==
%   max(X). Y must be an array of floating-point values whose size along the
%   first dimension is equal to the length of X, or a cell array of such arrays
%   such that all of Y{:} match the length of X. X_INTERVALS must be an Nx2
%   matrix of the same floating-point class as X, whose elements all fall within
%   the domain [X(1), X(end)]. Each 2-element row X_INTERVALS(i,:) specifies the
%   start and end of the ith interval.
%
%   X_TOL must be a positive scalar of the same floating-point class as X, which
%   specifies the desired density of resampling.
%
%   The return values X_RESAMPLED and Y_RESAMPLED are column cell arrays whose
%   length equals the number of time intervals in X_INTERVALS. Each cell of
%   X_RESAMPLED is a column vector of unique monotonically increasing values in
%   the domain, with no first-order differences greater than X_TOL, whose first
%   and last elements match the corresponding row of X_INTERVALS. Each cell of
%   Y_RESAMPLED is either an array or a nested cell array (depending on the
%   format of Y) which contains resampled function values.
%
%   RESAMPLE_ON_INTERVALS(...,METHOD) specifies alternate interpolation methods
%   ('linear' by default). Documentation for the MATLAB function INTERP1
%   describes alternative methods.
%
%   See also FIND_INTERVALS (written by smk).
%
%Written by smk 2009 June 24.
%

if isempty(x) || ~isvector(x) || ~isfloat(x) || ~isreal(x) || ...
    ~all(isfinite(x)) || ~all(diff(x) > 0) || ...
    (x(1) ~= min(x)) || (x(end) ~= max(x))
  error(['X must be a non-empty vector of real finite strictly ' ...
      'monotonically increasing floating-point values']);
end

if ~iscell(y)
  % flag for whether the output should be structured as an array
  y_cell_array_flag = false;
  y = {y};
else
  y_cell_array_flag = true;
end
if ~iscell(y) || ~all(cellfun(@(c) isfloat(c) && size(c,1)==numel(x),y))
  error(['Y must be a floating-point array whose length along the first ' ...
      'dimension equals the number of elements in X, or a cell array of ' ...
      'such arrays']);
end

if ~isnumeric(x_intervals) || (ndims(x_intervals) ~= 2) || ...
    ~isreal(x_intervals) || ~strcmp(class(x),class(x_intervals)) || ...
    any(x_intervals(:) < x(1)) || any(x_intervals(:) > x(end)) || ...
    ~all(x_intervals(:,1) < x_intervals(:,2))
  error(['X_INTERVALS must be an Nx2 matrix of the same floating-point ' ...
      'class as X, such that the values in the second column are all ' ...
      'greater than the values in the first column']);
end

if ~isnumeric(x_tol) || ~isreal(x_tol) || ~isscalar(x_tol) || ...
    ~strcmp(class(x),class(x_tol)) || ~(x_tol > 0)
  error('X_TOL must be a real positive scalar of the same class as X');
end

n_intervals = size(x_intervals,1);
x_resampled = cell([n_intervals 1]);
y_resampled = cell([n_intervals 1]);

for i = 1:n_intervals
  n_resample = 1 + ceil(diff(x_intervals(i,:))/x_tol);
  % apply UNIQUE to remove duplicates if x_tol small
  x_resampled{i} = unique( ...
      linspace(x_intervals(i,1),x_intervals(i,end),n_resample)');
  try
    if (y_cell_array_flag)
      y_resampled{i} = cellfun(@(arg) interp1(x,arg,x_resampled{i}, ...
          varargin{:}),y,'UniformOutput',false);
    else
      y_resampled{i} = interp1(x,y{1},x_resampled{i},varargin{:});
    end
  catch
    error('interpolation of data failed');
  end
end


