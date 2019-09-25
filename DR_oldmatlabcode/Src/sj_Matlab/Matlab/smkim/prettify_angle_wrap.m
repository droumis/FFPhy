function [x_pretty, a_pretty] = prettify_angle_wrap(x,a,center,tol)
%PRETTIFY_ANGLE_WRAP Insert NaN values and linearly interpolate data around phase-wrapping discontinuities for pretty plotting.
%
%   [X_PRETTY, A_PRETTY] = PRETTIFY_ANGLE_WRAP(X,A,CENTER,TOL), for column
%   vectors X and A of equal size, maps the elements of A (interpreted as
%   radians) to the interval [CENTER-pi, CENTER+pi], finds phase-wrapping jumps
%   greater than or equal to TOL between consecutive elements, inserts NaN
%   values into both X and A at these discontinuities, and linearly interpolates
%   before and after each discontinuity to achieve visually gapless coverage on
%   either side of the discontinuity.
%
%   PRETTIFY_ANGLE_WRAP(X,A,CENTER) assumes the default value TOL = pi.
%
%   See also INTERP1_ANGLE (written by smk), UNWRAP (MATLAB standard m-file).
%
%Depends on:
%   RECENTER_ANGLE (written by SMK)
%
%Written by smk, 2009 November 14.
%

if (exist('recenter_angle') ~= 2)
  error(['PRETTIFY_ANGLE_WRAP depends on m-file RECENTER_ANGLE ' ...
      '(written by smk)']);
end

% input checking
if ~isnumeric(x) || isempty(x) || (ndims(x) > 2) || (size(x,2) > 1) || ...
    ~isfloat(x) || ~isreal(x) || ~all(isfinite(x(:)))
  error(['X must be a non-empty column vector of real floating-point ' ...
      'finite values']);
end
if ~isnumeric(a) || isempty(a) || (ndims(a) > 2) || (size(a,2) > 1) || ...
    ~isfloat(a) || ~isreal(a)
  error('A must be a floating-point real column vector');
end
if any((a(:) < -pi) | (a(:) > pi))
  %warning('values of A outside the interval [-pi,+pi] will be wrapped');
end
if any(isnan(a(:))) || any(~isfinite(a(:)))
  warning('interpolation near NaN or Inf elements of A will be undefined');
end
% Check size compatibility between X and A.
if (numel(x) ~= numel(a))
  error('A and X must have the same number of elements');
end
if (nargin < 3)
  center = 0;
end
if (nargin < 4)
  tol = pi;
end
if ~isscalar(center) || ~isreal(center) || ~isfinite(center) || ...
    ~isfloat(center)
  error('CENTER must be a real finite floating-point scalar');
end
if ~isscalar(tol) || ~isreal(tol) || ~isfinite(tol) || ~isfloat(tol)
  error('TOL must be a real finite floating-point scalar');
end

% Convert to double, subtract offset, map to the interval [-pi,+pi], and unwrap
a_recentered = recenter_angle(double(a) - center,0);
try
  a_unwrapped = unwrap(a_recentered,tol);
catch
  error('TOL is not valid');
end
% Find jumps at -pi/+pi
jump_idx = find(abs(diff(a_recentered) - diff(a_unwrapped)) > tol);

% If no jumps, return original X and recentered A
if isempty(jump_idx)
  x_pretty = x;
  a_pretty = recenter_angle(double(a),center);
  a_pretty = cast(a_pretty,class(a));
  return;
end

% Otherwise, split the data at the jumps
n_intervals = numel(jump_idx)+1;
start_idx = [1; jump_idx+1];
end_idx = [jump_idx; numel(a)];
x_pretty = cell([n_intervals 1]);
a_pretty = cell([n_intervals 1]);

% Each cell contains a vector of values with a trailing NaN value
for i = 1:n_intervals
  j = start_idx(i):end_idx(i);
  if (i == 1)
    x_pretty{i} = [ NaN; x(j); NaN; NaN ];
    a_pretty{i} = center + ...
        [ pi*sign(a_recentered(j(1))); a_recentered(j); ...
        pi*sign(a_recentered(j(end))); NaN ];
    x_pretty{i}(end-1) = interp1( ...
        [ a_pretty{i}(end-2); ...
        a_pretty{i}(end-2) + diff(a_unwrapped([j(end) start_idx(i+1)])) ], ...
        [ x_pretty{i}(end-2); x(start_idx(i+1)) ], ...
        a_pretty{i}(end-1));
  elseif (i == n_intervals)
    x_pretty{i} = [ NaN; x(j); NaN];
    a_pretty{i} = center + ...
        [ pi*sign(a_recentered(j(1))); a_recentered(j); ...
        pi*sign(a_recentered(j(end))) ];
    x_pretty{i}(1) = interp1( ...
        [ a_pretty{i}(2) - diff(a_unwrapped([end_idx(i-1) j(1)])); ...
        a_pretty{i}(2) ], ...
        [ x(end_idx(i-1)); x_pretty{i}(2) ], ...
        a_pretty{i}(1));
  else
    x_pretty{i} = [ NaN; x(j); NaN; NaN];
    a_pretty{i} = center + ...
      [ pi*sign(a_recentered(j(1))); a_recentered(j); ...
        pi*sign(a_recentered(j(end))); NaN ];
    x_pretty{i}(1) = interp1( ...
        [ a_pretty{i}(2) - diff(a_unwrapped([end_idx(i-1) j(1)])); ...
        a_pretty{i}(2) ], ...
        [ x(end_idx(i-1)); x_pretty{i}(2) ], ...
        a_pretty{i}(1));
    x_pretty{i}(end-1) = interp1( ...
        [ a_pretty{i}(end-2); ...
        a_pretty{i}(end-2) + diff(a_unwrapped([j(end) start_idx(i+1)])) ], ...
        [ x_pretty{i}(end-2); x(start_idx(i+1)) ], ...
        a_pretty{i}(end-1));
  end
end

% Collapse cell arrays into long vectors
x_pretty = cell2mat(x_pretty);
a_pretty = cell2mat(a_pretty);

x_pretty = cast(x_pretty,class(x));
a_pretty = cast(a_pretty,class(a));

