function x_intervals = find_intervals(x,y,func,x_tol,varargin)
%FIND_INTERVALS Find intervals of domain on which a function satisfies a test condition
%   X_INTERVALS = FIND_INTERVALS(X,Y,FUNC,X_TOL) returns the set of
%   non-overlapping intervals on which the function Y, sampled at points X in
%   its real domain, satisfies the test condition specified by the boolean
%   function FUNC. X must be a column vector of strictly monotonically
%   increasing unique real numeric values, such that X(1) == min(X) and
%   X(end) == max(X). Y must be an array of floating-point values whose size
%   along the first dimension is equal to the length of X, or a cell array of
%   such arrays such that all of Y{:} match the length of X. FUNC must be a
%   function handle such that FUNC(Y{:}) returns a logical column vector of the
%   same size as X. If Y is a cell array, is up to the user to ensure that FUNC
%   accepts multiple arguments in the same order as they appear in the cell
%   array Y. X_INTERVALS is an Nx2 matrix in which each row is a real cadlag
%   interval [start, end).
%   
%   FIND_INTERVALS approximates the intervals of the domain [X(1) X(end)] on
%   which FUNC(Y) evaluates to true. It does this by finding consecutive ith,
%   (i+1)th samples such that FUNC(Y(i) ~= FUNC(Y(i+1)) and interpolating Y
%   between these samples on a grid whose spacing is less than or equal to X_TOL.
%   X is promoted to double-precision floating-point before finding intervals to
%   minimize round-off error. If Y is a cell array, the array in each cell is
%   interpolated independently. An error is raised if FUNC(Y) changes value more
%   than once between Y(i) and Y(i+1), so this function is not well-suited to
%   data in which FUNC(Y) fluctuates more frequently than the sampling density
%   of Y.
%
%   FIND_INTERVALS(...,METHOD) specifies alternate interpolation methods
%   ('linear' by default). Documentation for the MATLAB function INTERP1
%   describes alternative methods.
%
%   See also: IS_INTERVALS, RESAMPLE_ON_INTERVALS.
%
%Depends on:
%   IS_INTERVALS (written by SMK)
%
%Written by smk 2009 June 24.
%

if (exist('is_intervals') ~= 2)
  error(['FIND_INTERVALS depends on m-file IS_INTERVALS ' ...
      '(written by smk)']);
end

if isempty(x) || ~isvector(x) || ~isnumeric(x) || ~isreal(x) || ...
    ~all(isfinite(x)) || ~all(diff(x) > 0) || ...
    (x(1) ~= min(x)) || (x(end) ~= max(x))
  error(['X must be a non-empty vector of real finite strictly ' ...
      'monotonically increasing numeric values']);
end
% Convert X to double (will reconvert for output)
xclass = class(x);
x = double(x);

if ~iscell(y)
  y = {y};
end
if ~iscell(y) || ~all(cellfun(@(c) isfloat(c) && size(c,1)==numel(x),y))
  error(['Y must be a floating-point array whose length along the first ' ...
      'dimension equals the number of elements in X, or a cell array of ' ...
      'such arrays']);
end

if ~isa(func,'function_handle')
  error(['FUNC must be a function handle']);
end
try
  tf = func(y{:});
  assert(isa(tf,'logical') && isequal(size(tf),size(x)));
catch
  error(['FUNC must be a function that returns a logical vector of the ' ...
    'same dimensions as X']);
end

if ~isnumeric(x_tol) || ~isreal(x_tol) || ~isscalar(x_tol) || ...
    ~strcmp(xclass,class(x_tol)) || ~(x_tol > 0)
  error('X_TOL must be a real positive scalar of the same class as X');
end

switch nnz(tf)
case numel(x)
  x_intervals = zeros([1 2],xclass);
  x_intervals(1) = x(1);
  x_intervals(end) = x(end);
  return;
case 0
  warning('no intervals satisfy the test condition');
  x_intervals = zeros([0 2],xclass);
  return;
end

% Construct X grid in neighborhoods where FUNC(Y{:}) changes value between
% consecutive samples, and interpolate each cell of Y at points on this grid.
% The grid spacing is less than or equal to X_TOL.
xi = cell2mat(arrayfun(@(i) ...
    linspace(x(i),x(i+1),1 + ceil((x(i+1)-x(i))/x_tol))', ...
    find(diff(tf)),'UniformOutput',false));
try
  yi = cellfun(@(arg) interp1(x,arg,xi,varargin{:}),y,'UniformOutput',false);
catch
  error('error interpolating data');
end
try
  tfi = func(yi{:});
  assert(islogical(tfi) && isequal(size(tfi),size(xi)));
catch
  error('error evaluating FUNC on interpolated data');
end

% Find change-points
di = diff(tfi);
% Assemble column vectors of interval startpoints and endpoints
if tfi(1)
  x_start = [x(1); xi(1 + find(di > 0))];
else
  x_start = xi(1 + find(di > 0));
end
if tfi(end)
  x_end = [xi(find(di < 0)); x(end)];
else
  x_end = xi(find(di < 0));
end 
if ~isequal(x_start,x_end) && ~all(x_start <= x_end)
  error('there is a bug in FIND_INTERVALS. oops');
end
assert(isequal(size(x_start),size(x_end)));
% Construct output of desired class
x_intervals = zeros([numel(x_start) 2],xclass);
x_intervals(:,1) = x_start;
x_intervals(:,2) = x_end;

% Remove degenerate intervals of zero length
x_intervals(find(x_start == x_end),:) = [];

if isempty(x_intervals)
  warning('No intervals were found; returning empty intervals matrix');
end

if ~is_intervals(x_intervals)
  error('There is a bug in either FIND_INTERVALS or IS_INTERVALS');
end


