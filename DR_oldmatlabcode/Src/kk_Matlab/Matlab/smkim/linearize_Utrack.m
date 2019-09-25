function linpos = linearize_Utrack(rawpos,transform,params)
%
% This function performs an automated linearization and smoothing given raw
% position data. THERE ARE HARD-CODED CONSTANTS IN THIS FUNCTION!!!!
%
% The output is a linpos struct, which inherits all of the input pos struct and
% adds more fields
% 
%
% typical values for linearized smoothing params:
%   front_back_marker_weights: [1 0]
%   loess_halfwidth: 0.75
%   num_loess_iterations: 3
%   rrm_halfwidths: [0.1 0.1 0.1]
%   roughness_tolerance: 0.1 (cm)
%   epsilon: 2
%   min_stop_duration: 0.5
%
%Depends on:
%   TFORMFWD (MATLAB Image Processing Toolbox)
%   TFORMWD (MATLAB Image Processing Toolbox)
%   PDIST (MATLAB Statistics Toolbox)
%   FIND_NEARBY_MEX (written by smk)
%   IS_RAWPOS (written by smk)
%   IS_LINEARIZED (written by smk)
%   BISQUARE (written by smk)
%   TRICUBE (writtten by smk)
%   SPAPS (MATLAB Spline Toolbox)
%
%Written by smk (2009 November 1)
%
%

if (exist('tformfwd') ~= 2)
  error(['LINEARIZE_UTRACK depends on m-file TFORMFWD ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('tforminv') ~= 2)
  error(['LINEARIZE_UTRACK depends on m-file TFORMINV ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('findbounds') ~= 2)
  error(['LINEARIZE_UTRACK depends on m-file FINDBOUNDS ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('pdist') ~= 2)
  error(['LINEARIZE_UTRACK depends on m-file PDIST ' ...
      '(MATLAB Statistics Toolbox)']);
end
if (exist('find_nearby_mex') ~= 3)
  error(['LINEARIZE_UTRACK depends on mex-file FIND_NEARBY_MEX ' ...
      '(written by smk)']);
end
if (exist('is_rawpos') ~= 2)
  error('LINEARIZE_UTRACK depends on m-file IS_RAWPOS (written by smk)');
end
if (exist('is_linearized') ~= 2)
  error('LINEARIZE_UTRACK depends on m-file IS_LINEARIZED (written by smk)');
end
if (exist('bisquare') ~= 2)
  error('LINEARIZE_UTRACK depends on m-file BISQUARE (written by smk)');
end
if (exist('tricube') ~= 2)
  error('LINEARIZE_UTRACK depends on m-file TRICUBE (written by smk)');
end
if (exist('spaps') ~= 2)
  error(['LINEARIZE_UTRACK depends on m-file SPAPS ' ...
      '(MATLAB Spline Toolbox)']);
end

% DO NOT TAMPER WITH THIS CONSTANT CONVERSION FACTOR!
TS_PER_SEC = 1e4;

if ~is_rawpos(rawpos)
  error('RAWPOS does not appear to contain valid video head-tracking data');
end
% Check that transform struct is well-formed and has the same number of elements
% as rawpos
if ~isstruct(transform) || ~all(isfield(transform, ...
    {'ndims_in','ndims_out','forward_fcn','inverse_fcn','tdata'})) || ...
    ~isequal(transform.ndims_in,transform.ndims_out,2)
  error(['TRANSFORM does not appear to be a valid image transform struct ' ...
      'of the type returned by CP2TFORM']);
end
if (numel(transform) ~= numel(rawpos))
  error(['number of struct elements of TRANSFORM does not match number ' ...
      'of struct elements of RAWPOS']);
end
% Check that params struct is well-formed
if ~isstruct(params) || (numel(params) ~= 1)
  error('smoothing params must be a single struct');
end
REQUIRED_PARAMS = {
    'front_back_marker_weights'    , ...
    'loess_halfwidth'              , ...
    'num_loess_iterations'         , ...
    'rrm_halfwidths'               , ...
    'min_stop_duration'            , ...
    'roughness_tolerance'            , ...
    'epsilon'                      , ...
    'min_front_back_separation'        };
if ~all(isfield(params,REQUIRED_PARAMS))
  error('params is missing one or more required fields');
end

param_names = fieldnames(params);
for i = 1:numel(param_names)
  name = param_names{i};
  value = params.(name);
  switch name
  case 'front_back_marker_weights'
    if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
        ~(numel(value) == 2) || ~all(value >= 0) || (sum(value) ~= 1)
      error('%s must be a 2-element vector of non-negative reals',name);
    end
  case 'loess_halfwidth'
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
      ~isfinite(value) || (value <= 0)
      error('%s must be a real positive scalar',name);
    end
  case 'num_loess_iterations'
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || (value < 1) || (round(value) ~= value)
      error('%s must be a real positive integer',name);
    end
  case 'rrm_halfwidths'
    if ~isnumeric(value) || ~isreal(value) || ~isvector(value) || ...
        (numel(value) < 2) || ~all(isfinite(value)) || ~all(value > 0) || ...
        ~all(diff(value) <= 0)
      error('%s must be a real vector of decreasing positive values',name);
    end
  case 'min_stop_duration'
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || (value <= 0)
      error('%s must be a real positive scalar',name);
    end
  case 'roughness_tolerance'
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || (value <= 0)
      error('%s must be a positive real float scalar',name);
    end
  case 'epsilon'
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || (value <= 0)
      error('%s must be a positive real float scalar',name);
    end
  case 'min_front_back_separation'
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        ~isfinite(value) || (value <= 0)
      error('%s must be a real non-negative scalar',name);
    end
  case 'head_direction_epsilon'
    if ~isnumeric(value) || ~isreal(value) || ~isscalar(value) || ...
        (value < 0)
      error('%s must be a real positive scalar',name);
    end
  end
end

% Initialize POSITION and inherit meta-data from RAWPOS
linpos = rmfield(rawpos,{'timestamp','xfront','yfront','xback','yback'});
[linpos(:).units] = deal('centimeters');
[linpos(:).params] = deal(params);

for i = 1:numel(rawpos)
  % Determine whether a back marker was tracked
  if isempty(rawpos(i).xback) && isempty(rawpos(i).yback)
    back_exists = false;
  elseif ~isempty(rawpos(i).xback) && ~isempty(rawpos(i).yback)
    back_exists = true;
  else
    error('xback and yback must be the same size');
  end  
  % Convert timestamps to floating-point seconds
  t = double(rawpos(i).timestamp)/TS_PER_SEC;
  % Censor samples with NaN or zero values
  include_mask = ...
      ~isnan(rawpos(i).xfront) & (rawpos(i).xfront ~= 0) & ...
      ~isnan(rawpos(i).yfront) & (rawpos(i).yfront ~= 0);
  if (back_exists)
    include_mask = include_mask & ...
        ~isnan(rawpos(i).xback) & (rawpos(i).xback ~= 0) & ...
        ~isnan(rawpos(i).yback) & (rawpos(i).yback ~= 0);
  end
  % Censor and promote to double to take advantage of MATLAB acceleration
  t_censored = t(include_mask);
  ufront_censored = double(rawpos(i).xfront(include_mask));
  vfront_censored = double(rawpos(i).yfront(include_mask));
  if (back_exists)
    uback_censored = double(rawpos(i).xback(include_mask));
    vback_censored = double(rawpos(i).yback(include_mask));
  end
  % Transform from image coordinates to real-world coordinates
  if isa(transform(i).forward_fcn,'function_handle')
    [xfront_censored, yfront_censored] = ...
        tformfwd(transform(i),ufront_censored,vfront_censored);
    if (back_exists)
      [xback_censored, yback_censored] = ...
          tformfwd(transform(i),uback_censored,vback_censored);
    end
  else
    % construct lookup table in the inverse direction, then interpolate in the
    % forward direction
    if (back_exists)
      bounds = findbounds(transform(i), ...
          [min(min(ufront_censored,uback_censored)), ...
          min(min(vfront_censored,vback_censored)); ...
          max(max(ufront_censored,uback_censored)), ...
          max(max(vfront_censored,vback_censored))]);
    else
      bounds = findbounds(transform(i), ...
          [min(ufront_censored), min(vfront_censored); ...
          max(ufront_censored), max(vfront_censored)]);
    end
    [xgrid, ygrid] = meshgrid( ...
        floor(bounds(1,1) - 2):ceil(bounds(2,1) + 2), ...
        floor(bounds(1,2) - 2):ceil(bounds(2,2) + 2) );
    [ugrid, vgrid] = tforminv(transform(i),xgrid,ygrid);
    % x and y are expressed in real-world coordinates with units of centimeters
    xfront_censored = griddata(ugrid,vgrid,xgrid, ...
        ufront_censored,vfront_censored,'cubic');
    yfront_censored = griddata(ugrid,vgrid,ygrid, ...
        ufront_censored,vfront_censored,'cubic');
    if (back_exists)
      xback_censored = griddata(ugrid,vgrid,xgrid, ...
          uback_censored,vback_censored,'cubic');
      yback_censored = griddata(ugrid,vgrid,ygrid, ...
          uback_censored,vback_censored,'cubic');
    end
  end
  if (back_exists)
    % Estimate head position from weighted average of front and back
    if (size(params.front_back_marker_weights,1) > 1)
      % Coerce to be a 2-element row vector for proper alignment
      params.front_back_marker_weights = params.front_back_marker_weights';
    end
    x_censored = sum(bsxfun(@times,[xfront_censored xback_censored], ...
        params.front_back_marker_weights),2) ./ ...
        sum(params.front_back_marker_weights);
    y_censored = sum(bsxfun(@times,[yfront_censored yback_censored], ...
        params.front_back_marker_weights),2) ./ ...
        sum(params.front_back_marker_weights);
  else
    x_censored = xfront_censored;
    y_censored = yfront_censored;
  end
  
  % Keep a copy of the transform
  linpos(i).transform = transform(i);
  % Map (x,y) coordinates to linearized position along track
  lin_censored = cart2lin(x_censored,y_censored);

  % Perform robust iterative LOESS smoothing to get outlier-resistant denoised
  % linearized position and weights
  loess = loess_trajectory(t_censored,lin_censored,t, ...
      params.loess_halfwidth,params.num_loess_iterations);
  % Interpolate LOESS fit to estimate first and second derivatives and use these
  % to locally adjust smoothing spline
  sp_interp = spaps(t,loess,0,3);
  D1 = fnval(fnder(sp_interp,1),t);
  D2 = fnval(fnder(sp_interp,2),t);
  assert(~any(isnan(D1)) && ~any(isnan(D2)));
  ti = 0.5*(t(2:end) + t(1:end-1));
  i_nearby = find_nearby_mex(ti,t, ...
      -params.loess_halfwidth,+params.loess_halfwidth);
  % TOL argument of SPAPS is vector such that first element TOL(1) specifies the
  % constraint on the total sum of the squared residuals over all points, and
  % TOL(2:end) specify the adaptive weighting of the roughness penalty on the
  % intervals between the points.
  tol = zeros([numel(loess), 1]);
  tol(1) = numel(t) * params.roughness_tolerance^2;
  for j = 1:numel(ti)
    lags = t(i_nearby{j}) - ti(j); 
    weights = tricube(lags/params.loess_halfwidth);
    % Ratio of (kernel-weighted integrated squared second derivative) to
    % (kernel-weighted integrated squared first derivative). If the LOESS
    % estimate is very wiggly, then this ratio will be large. Ratios smaller
    % than unity are truncated at unity (otherwise we end up overfitting).
    tol(j+1) = max(1,sum(weights .* D2(i_nearby{j}).^2) / ...
        sum(weights .* D1(i_nearby{j}).^2));
  end
  % Apply smoothing-spline to the LOESS estimate with a locally-adapted
  % roughness penalty. We compute the roughness penalty on the 3rd derivative in
  % order to obtain smooth first derivatives. SPAPS ignores NaN values.
  sp = spaps(t,loess,tol,3);
  smooth_position = fnval(sp,t);
  smooth_speed = fnval(fnder(sp,1),t);

  % Reverse mapping to get idealized (x,y) coordinates projected onto the track
  [x_projection, y_projection, track_direction] = lin2cart(smooth_position);
  % Find arrests by applying repeated running median filter to both the loess
  % smooth and the final smooth
  arrests = union_intervals( ...
      find_intervals(uint32(TS_PER_SEC*t), ...
      double(rrm_arrests(t,smooth_position, ...
      params.rrm_halfwidths,params.min_stop_duration,params.epsilon)), ...
      @(x) x > 0.5,uint32(1)), ...
      find_intervals(uint32(TS_PER_SEC*t(~isnan(loess))), ...
      double(rrm_arrests(t(~isnan(loess)),loess(~isnan(loess)), ...
      params.rrm_halfwidths,params.min_stop_duration,params.epsilon)), ...
      @(x) x > 0.5,uint32(1)) );
  % Join together arrests that are separated by less than
  % 2*max(params.rrm_halfwidths)
  arrests = intersect_intervals(rawpos(i).timestamp([1 end])', ...
      join_intervals(arrests,uint32(TS_PER_SEC*2*max(params.rrm_halfwidths))));

  % Assign to output
  linpos(i).timestamp = rawpos(i).timestamp;
  linpos(i).x_projection = single(x_projection);
  linpos(i).y_projection = single(y_projection);
  linpos(i).position = single(smooth_position);
  linpos(i).speed = single(smooth_speed);
  linpos(i).arrests = arrests;
end

if ~is_linearized(linpos)
  error('There is bug in either LINEARIZE_UTRACK or IS_LINEARIZED');
end

end % end main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOESS ITERATIVE NONLINEAR SMOOTHER TO ESTIMATE LINEARIZED POSITION AND
% VELOCITY
function [smooth_val, smooth_der] = ...
    loess_trajectory(t,lin,ti,halfwidth,num_iterations)
  if any(~isfinite(lin)) || any(~isfinite(t)) || any(diff(t) <= 0)
    error(['input vectors must be finite non-NaN, and t vector must be ' ...
        'monotonically increasing']);
  end
  % Degree of the local smoothing polynominal: linear is bad, quadratic is
  % adequate, cubic tends to overfit.
  LOESS_DEGREE = 2;
  n = numel(ti);
  smooth_val = nan([n 1]);
  smooth_der = nan([n 1]);
  % Turn off warning in the LSCOV function, as otherwise MATLAB spits out a
  % flood of warning statements in the console
  warning off MATLAB:lscov:RankDefDesignMat
  % Index look-up for nearby samples that fall within the smoothing bandwidth
  i_nearby = find_nearby_mex(ti,t,-halfwidth,+halfwidth);
  % Keep track of points where the loess fit fails because of missing samples
  % in the neighborhood
  failpoints = false([n 1]);
  % Perform iterative loess at each position sample
  coeffs = zeros([n LOESS_DEGREE+1]);
  for i = 1:n
    % Construct degree-d Vandermonde matrix: first column is dth-order
    % coefficient, second column is (d-1)th order coefficient, ... dth column is
    % linear term coefficient, (d+1)th column is constant term (this is the
    % convention for POLYVAL, POLYDER, POLYFIT, etc.)
    lags = t(i_nearby{i}) - ti(i); 
    basis = ones([numel(lags) LOESS_DEGREE+1]);
    for d = LOESS_DEGREE:-1:1
      basis(:,d) = lags .* basis(:,d+1);
    end
    % Grab observations in this window
    obs = lin(i_nearby{i});
    % Initialize regression weights by tricube weighting over the smoothing
    % bandwidth
    weights = tricube(lags/halfwidth);
    iter_count = 0;
    while true
      % Fit polynomial to local observations (weighted by weights vector)
      try
        coeffs(i,:) = lscov(basis,obs,weights);
      catch
        % LSCOV will fail if there are not enough neighboring data points
        % within the bandwidth, in which case we flag this sample and move on
        % to the next sample
        failpoints(i) = true;
        break;
      end
      if (iter_count < num_iterations)
        % Adjust the weights acccording to the residuals: normalize residuals
        % w.r.t. 6*median(residuals), then update weights by multiplying with
        % bi-square function of normalized residuals
        residuals = polyval(coeffs(i,:),lags) - obs;
        weights = weights .* ...
            bisquare(residuals/(6*median(abs(residuals))));
        iter_count = iter_count + 1;
      else
        break;
      end
    end
  end
  % Now that the smoothing is complete, estimate linearized position from the
  % final fit coefficients. We leave NaN values where the smooth could not be
  % determined
  smooth_val(~failpoints) = coeffs(find(~failpoints),end); 
  smooth_der(~failpoints) = coeffs(find(~failpoints),end-1); 
end % end function LOESS_TRAJECTORY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION FOR MAPPING (X,Y) COORDINATES TO LINEARIZED POSITIONS
% Hard-coded constants to match dimensions of my U track.
function lin = cart2lin(x,y);
  lin = nan(size(x));

  domain{1} = (y >= 22) & (x < 22);
  lin(domain{1}) = -22 - y(domain{1});

  domain{2} = (y < 22) & (y >= x);
  lin(domain{2}) = 22*((22 - y(domain{2}))./(22 - x(domain{2})) - 2);

  domain{3} = (y < 22) & (y < x) & (y < 44 - x);
  lin(domain{3}) = 22*(x(domain{3}) - 22)./(22 - y(domain{3}));

  domain{4} = (y < 22) & (y >= 44 - x);
  lin(domain{4}) = 22*(2 - (22 - y(domain{4}))./(x(domain{4}) - 22));

  domain{5} = (y >= 22) & (x >= 22);
  lin(domain{5}) = 22 + y(domain{5});

  if any(xor(isnan(lin),(isnan(x) | isnan(y))))
    error('there is a mistake in cart2lin');
  end
end % end function CART2LIN


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION FOR MAPPING LINEARIZED POSITIONS TO IDEALIZED (X,Y) COORDINATES
% Hard-coded constants to match dimensions of my U track.
function [x, y, track_direction] = lin2cart(lin)
  x = nan(size(lin));
  y = nan(size(lin));
  % track_direction is the idealized direction in which linearized position
  % increases along the track
  track_direction = nan(size(lin));

  domain{1} = (lin <= -22);
  x(domain{1}) = 0;
  y(domain{1}) = -22 - lin(domain{1});
  track_direction(domain{1}) = -pi/2;

  domain{2} = (lin > -22) & (lin < 22);
  x(domain{2}) = lin(domain{2});
  y(domain{2}) = 0;
  track_direction(domain{2}) = 0;

  domain{3} = (lin >= 22);
  x(domain{3}) = 44;
  y(domain{3}) = lin(domain{3}) - 22;
  track_direction(domain{3}) = +pi/2;
  
  if any(xor(isnan(lin),(isnan(x) | isnan(y))))
    error('there is a mistake in lin2cart');
  end
end % end function LIN2CART


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPEATED RUNNING MEDIAN DETECTION OF STOPS
function [stopped, lin_median] = ...
    rrm_arrests(t,lin,halfwidths,min_stop_duration,epsilon)
  n = numel(t);
  % Sequentially compute running medians over x and y, iterating through the
  % halfwidths specified in params.rrm_halfwidths
  lin_median = nan(size(lin));
  for hw = halfwidths
    % i_nearby{i} contains indices to all samples that are within halfwidth of
    % the ith sample
    i_nearby = find_nearby_mex(t,t,-hw,+hw);
    % Surprisingly, this for loop is much faster than arrayfun
    for i = 1:n
      lin_median(i) = median(lin(i_nearby{i}));
    end
    % Replace lin with its running median
    lin = lin_median;
  end
  % After applying the running median smoothers, identify continuous intervals
  % that satisfy the criteria min_stop_duration (proximity in time) and epsilon
  % (proximity in space). For each sample, look ahead to min_stop_duration later
  % and determine whether all samples in this range fall within epsilon
  stopped = false([n 1]);
  for i = 1:(n-1)
    i_range = i : ...
        (i + find(t((i+1):end) - t(i) >= min_stop_duration,1,'first'));
    if ~isempty(i_range)
      stopped(i_range) = stopped(i_range) | ...
          all(pdist(lin(i_range)) <= epsilon);
    end
  end
end % end function RRM_ARRESTS

%{
% For documentation of hard-coded U-track dimensions:

% SEGMENTS is a struct array in which each element corresponds to a linear
% segment the Utrack. The user selects points along each edge, which are then
% interpolated/resampled to produce control points for image registration. The
% start->end order of traversal along each edge is important. The fields xstart,
% ystart, xend, yend are expressed in centimeters relative to the location of
% one of the food wells.
SEGMENTS = cell2struct({ ...
    0   , 147 , 0   , 0 ; ...
    0   , 0   , 44  , 0 ; ...
    44  , 0   , 44  , 147   }, ...
    {'xstart','ystart','xend','yend'},2);

% struct array of linearization functions
linearization_mapping(1) = struct( ...
    'select_fcn'  , @(x,y) (y >= 22) & (x < 22)               , ...
    'x_fcn'       , @(x,y) zeros(size(x),class(x))            , ...
    'y_fcn'       , @(x,y) y                                  , ...
    'dist_fcn'    , @(x,y) -22 - y                            );
linearization_mapping(2) = struct( ...
    'select_fcn'  , @(x,y) (y < 22) & (y >= x)                , ...
    'x_fcn'       , @(x,y) zeros(size(x),class(x))            , ...
    'y_fcn'       , @(x,y) 22*(1 - (22 - y)./(22 - x))        , ...
    'dist_fcn'    , @(x,y) 22*((22 - y)./(22 - x) - 2)        );
linearization_mapping(3) = struct( ...
    'select_fcn'  , @(x,y) (y < 22) & (y < x) & (y < 44 - x)  , ...
    'x_fcn'       , @(x,y) 22*(1 + .x - 22)./(22 - y))        , ...
    'y_fcn'       , @(x,y) zeros(size(y),class(y))            , ...
    'dist_fcn'    , @(x,y) 22*(x - 22)./(22 - y)              );
linearization_mapping(4) = struct( ...
    'select_fcn'  , @(x,y) (y < 22) & (y >= 44 - x)           , ...
    'x_fcn'       , @(x,y) 44*ones(size(s.x),class(s.x))      , ...
    'y_fcn'       , @(x,y) 22*(1 - (22 - y)./(x - 22))        , ...
    'dist_fcn'    , @(x,y) 22*(2 - (22 - y)./(x - 22))        );
linearization_mapping(5) = struct( ...
    'select_fcn'  , @(x,y) (y >= 22) & (x >= 22)              , ...
    'x_fcn'       , @(x,y) 44*ones(size(x),class(x))          , ...
    'y_fcn'       , @(x,y) y                                  , ...
    'dist_fcn'    , @(x,y) 22 + y                             );
%}

