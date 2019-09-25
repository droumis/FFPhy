function position = estimate_position(rawpos,transform,params)
%ESTIMATE_POSITION Estimate head direction, smoothed position and velocity from raw video pixel coordinates of head-tracking LEDs.
%
%   POSITION = ESTIMATE_POSITION(RAWPOS,TRANSFORM,SMOOTHING_PARAMS) estimates
%   head direction, smoothed position and velocity from the coordinates of the
%   head-tracking LEDs in the RAWPOS.
%
%   First, the video pixel coordinates in RAWPOS are transformed to real-world
%   allocentric coordinates via the transform specified by TRANSFORM. TRANSFORM
%   must be an array of structs of the type returned by the CP2TFORM in the
%   MATLAB Image Processing Toolbox. The transformation is applied by calling
%   TFORMFWD in the MATLAB Image Processing Toolbox. TRANSFORM must be of the
%   same size as RAWPOS; the pixel coordinates in RAWPOS(i) are transformed
%   accordingly to TRANSFORM(i).
%
%   Next, the line segment is drawn connecting the positions of front and back
%   LEDs in the real-world coordinate frame (e.g. units of centimeters). The
%   direction of this line segment is filtered to estimate head direction. A
%   point along this line segment (between the two LEDs) is chosen and its
%   location is smoothed with LOESS regression. Stopping times are identified
%   by applying a repeated running median smoother to the smoothed position to
%   identify stops (note that "stopped" does not necessarily mean immobile or
%   asleep; whisking, grooming, etc. may be occuring).
%
%   PARAMS is a struct that contains the following fields:
%
%   front_back_marker_weights (applicable only when both back and front markers
%   are tracked)
%     Parameter for the relative weighting of the front LEDs versus back LEDs
%     when estimating "true" head position. This is a two-element vector of
%     positive reals. Typical value is [1 0] (i.e. rely on the front LEDs).
%
%   loess_halfwidth*
%     The halfwidth (expressed in seconds) of the sliding window on which
%     LOESS is done. This is the most important parameter for estimating
%     smooth position and velocity. Typical value is 0.75 (750 ms).
%
%   num_loess_iterations*
%     Number of reweighting iterations to apply LOESS for robust rejection of
%     outliers. Unless the data are especially noisy, the estimate practically
%     converges after 3 or 4 iterations.
%
%   rrm_halfwidths*
%     Vector of halfwidths (in seconds) of sliding windows in which running
%     medians are successively computed for the repeated running median
%     smoother. Data are smoothed with a running median of halfwidth
%     rrm_halfwidths(1), and then this result is smoothed with a running median
%     of rrm_halfwidths(2), etc. Window sizes must be non-increasing:
%     rrm_halfwidths(1) >= rrm_halfwidths(2) >= ... >= rrm_halfwidths(end).
%     Typical value is [0.1 0.1 0.1].
%
%   epsilon*
%     The "closeness" threshold (in cm) for determining when the subject is
%     immobile. This should be a very small value, as small as or smaller than
%     the pixel resolution. Typical value is 0.5.
%
%   min_stop_duration*
%     Define the minimum duration (in seconds) of consecutive "close" position
%     samples that must occur for an interval to qualify as a stop. Typical
%     value is 0.5 (500 ms).
%
%   min_front_back_separation (applicable only when both back and front markers
%   are tracked)
%     When the pitch of the rat's head is excessive, the front and back lights
%     are not distinguishable any more; therefore, head direction estimates are
%     unreliable at these times. This (float) parameter is the minimum apparent
%     front-back separation (in cm) required for estimating head direction.
%     Obviously, the value depends on the size of the LED tracker boom.
%
%   *The smoothing parameters that are marked with asterisk are described in
%   Hen I., Sakov A., Kafkafi N., Golani I., Benjamini Y. (2004) The dynamics
%   of spatial behavior: how can robust smoothing techniques help? _Journal of
%   Neuroscience Methods_ 133: 161-172.
%
%   For explanation of regression model that is used to estimate head
%   direction, see section 6.4 of Fisher N.I. (1996)  _Statistical Analysis of
%   Circular Data_. Cambridge UP.
%
%Depends on:
%   TFORMFWD (MATLAB Image Processing Toolbox)
%   TFORMWD (MATLAB Image Processing Toolbox)
%   PDIST (MATLAB Statistics Toolbox)
%   FIND_NEARBY_MEX (written by smk)
%   IS_RAWPOS (written by smk)
%   IS_POSITION (written by smk)
%   BISQUARE (written by smk)
%   TRICUBE (writtten by smk)
%
%Written by smk (2009 November 1)
%

if (exist('tformfwd') ~= 2)
  error(['ESTIMATE_POSITION depends on m-file TFORMFWD ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('tforminv') ~= 2)
  error(['ESTIMATE_POSITION depends on m-file TFORMINV ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('findbounds') ~= 2)
  error(['ESTIMATE_POSITION depends on m-file FINDBOUNDS ' ...
      '(MATLAB Image Processing Toolbox)']);
end
if (exist('pdist') ~= 2)
  error(['ESTIMATE_POSITION depends on m-file PDIST ' ...
      '(MATLAB Statistics Toolbox)']);
end
if (exist('find_nearby_mex') ~= 3)
  error(['ESTIMATE_POSITION depends on mex-file FIND_NEARBY_MEX ' ...
      '(written by smk)']);
end
if (exist('is_rawpos') ~= 2)
  error('ESTIMATE_POSITION depends on m-file IS_RAWPOS (written by smk)');
end
if (exist('is_position') ~= 2)
  error('ESTIMATE_POSITION depends on m-file IS_POSITION (written by smk)');
end
if (exist('bisquare') ~= 2)
  error('ESTIMATE_POSITION depends on m-file BISQUARE (written by smk)');
end
if (exist('tricube') ~= 2)
  error('ESTIMATE_POSITION depends on m-file TRICUBE (written by smk)');
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
if ~isstruct(params) || ~isscalar(params)
  error('smoothing params must be a scalar struct');
end
REQUIRED_PARAMS = {
    'front_back_marker_weights'        , ...
    'loess_halfwidth'                 , ...
    'num_loess_iterations'            , ...
    'rrm_halfwidths'                   , ...
    'min_stop_duration'            , ...
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
position = rmfield(rawpos,{'timestamp','xfront','yfront','xback','yback'});
[position(:).units] = deal('centimeters');
[position(:).params] = deal(params);

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
  % Promote to double to take advantage of MATLAB acceleration
  ufront = double(rawpos(i).xfront);
  vfront = double(rawpos(i).yfront);
  uback = double(rawpos(i).xback);
  vback = double(rawpos(i).yback);
  % Censor samples with NaN or zero values
  censor_mask = ...
      isnan(ufront) | (ufront == 0) | isnan(vfront) | (vfront == 0);
  if (back_exists)
    censor_mask = censor_mask | ...
      isnan(uback) | (uback == 0) | isnan(vback) | (vback == 0);
  end
  t(censor_mask) = [];
  ufront(censor_mask) = [];
  vfront(censor_mask) = [];
  if (back_exists)
    uback(censor_mask) = [];
    vback(censor_mask) = [];
  end
  % n is the number of elements in t *after* NaN samples are removed
  n = numel(t);
  % Transform from image coordinates to real-world coordinates
  if isa(transform(i).forward_fcn,'function_handle')
    [xfront, yfront] = tformfwd(transform(i),ufront,vfront);
    if (back_exists)
      [xback, yback] = tformfwd(transform(i),xback,yback);
    end
  else
    % construct lookup table in the inverse direction, then interpolate in the
    % forward direction
    bounds = findbounds(transform(i), ...
        [min(min(ufront,uback)), min(min(vfront,vback)); ...
        max(max(ufront,uback)), max(max(vfront,vback))]);
    [xgrid, ygrid] = meshgrid( ...
        floor(bounds(1,1) - 2):ceil(bounds(2,1) + 2), ...
        floor(bounds(1,2) - 2):ceil(bounds(2,2) + 2) );
    [ugrid, vgrid] = tforminv(transform(i),xgrid,ygrid);
    % x and y are expressed in real-world coordinates with units of centimeters
    xfront = griddata(ugrid,vgrid,xgrid,ufront,vfront,'cubic');
    yfront = griddata(ugrid,vgrid,ygrid,ufront,vfront,'cubic');
    if (back_exists)
      xback = griddata(ugrid,vgrid,xgrid,uback,vback,'cubic');
      yback = griddata(ugrid,vgrid,ygrid,uback,vback,'cubic');
    end
  end
  % Keep a copy of the transform
  position(i).transform = transform(i);
  % Inherit timestamps
  position(i).timestamp = rawpos(i).timestamp(~censor_mask);
  if (back_exists)
    % Estimate head position from weighted average of front and back
    if (size(params.front_back_marker_weights,1) > 1)
      % Coerce to be 2-element row vector for proper alignment
      params.front_back_marker_weights = params.front_back_marker_weights';
    end
    x = sum(bsxfun(@times,[xfront xback], ...
        params.front_back_marker_weights),2) ./ ...
        sum(params.front_back_marker_weights);
    y = sum(bsxfun(@times,[yfront yback], ...
        params.front_back_marker_weights),2) ./ ...
        sum(params.front_back_marker_weights);
  else
    x = xfront;
    y = yfront;
  end

  % Perform LOESS smoothing to get denoised position and velocity
  [x_smooth, y_smooth, x_vel, y_vel] = loess_trajectory( ...
      t,x,y, ...
      params.loess_halfwidth,params.num_loess_iterations);
  position(i).x = single(x_smooth);
  position(i).y = single(y_smooth);
  position(i).x_velocity = single(x_vel);
  position(i).y_velocity = single(y_vel);
  position(i).x_residual = single(x - x_smooth);
  position(i).y_residual = single(y - y_smooth);
  % Perform LOESS smoothing to get denoised head direction
  if (back_exists)
    position(i).head_direction = single(loess_head_direction( ...
        t,xfront,yfront,xback,yback, ...
        params.loess_halfwidth,params.num_loess_iterations, ...
        params.min_front_back_separation));
  else
    position(i).head_direction = nan(size(t),'single');
  end  
  % Find arrests by applying repeated running median filter to the smoothed
  % trajectory
  position(i).stopped = rrm_arrests( ...
      t,x_smooth,y_smooth, ...
      params.rrm_halfwidths,params.min_stop_duration,params.epsilon);

  disp(position(i));
end

if ~is_position(position)
  error('There is bug in either ESTIMATE_POSITION or IS_POSITION');
end

end % end main function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOESS ITERATIVE NONLINEAR SMOOTHER TO ESTIMATE POSITION AND VELOCITY
function [x_smooth, y_smooth, x_vel, y_vel] = ...
      loess_trajectory(tuncensored,t,x,y,halfwidth,num_iterations)

  if any(~isfinite(x)) || any(~isfinite(y)) || any(~isfinite(t)) || ...
      any(diff(t) <= 0)
    error(['input vectors must be finite non-NaN, and t vector must be ' ...
        'monotonically increasing']);
  end
  % Degree of the local smoothing polynominal; quadratic is adequate
  LOESS_DEGREE = 2;
  n = numel(t);
  % Turn off warning in the LSCOV function, as otherwise MATLAB spits out a
  % flood of warning statements in the console
  warning off MATLAB:lscov:RankDefDesignMat
  % Index look-up for nearby samples that fall within the smoothing bandwidth
  i_nearby = find_nearby_mex(tuncensored,t,-halfwidth,+halfwidth);
  % Keep track of points where the loess fit fails because of missing samples
  % in the neighborhood
  failpoints = false([n 1]);
  % Perform iterative loess at each position sample
  x_coeffs = zeros([n LOESS_DEGREE+1]);
  y_coeffs = zeros([n LOESS_DEGREE+1]);
  for i = 1:n
    % Construct degree-d Vandermonde matrix: first column is dth-order
    % coefficient, second column is (d-1)th order coefficient, ... dth column is
    % linear term coefficient, (d+1)th column is constant term (this is the
    % convention for POLYVAL, POLYDER, POLYFIT, etc.)
    lags = t(i_nearby{i}) - t(i); 
    basis = ones([numel(lags) LOESS_DEGREE+1]);
    for d = LOESS_DEGREE:-1:1
      basis(:,d) = lags .* basis(:,d+1);
    end
    % Grab observations in this window
    x_obs = x(i_nearby{i});
    y_obs = y(i_nearby{i});
    % Initialize regression weights by tricube weighting over the smoothing
    % bandwidth
    weights = tricube(lags/halfwidth);
    iter_count = 0;
    while true
      % Fit polynomial to local observations (weighted by weights vector)
      try
        x_coeffs(i,:) = lscov(basis,x_obs,weights);
        y_coeffs(i,:) = lscov(basis,y_obs,weights);
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
        residuals = hypot(polyval(x_coeffs(i,:),lags) - x_obs, ...
            polyval(y_coeffs(i,:),lags) - y_obs);
        weights = weights .* ...
            bisquare(residuals/(6*median(abs(residuals))));
        iter_count = iter_count + 1;
      else
        break;
      end
    end
  end
  % Now that the smoothing is complete, estimate position and velocity from the
  % final fit coefficients
  x_smooth = x_coeffs(:,end); 
  y_smooth = y_coeffs(:,end); 
  x_vel = x_coeffs(:,end-1); 
  y_vel = y_coeffs(:,end-1);
  % Patch up the points where polynomial fit failed
  % Guarantee endpoint assignments for interpolation
  if failpoints(1)
    x_smooth(1) = x(1);
    y_smooth(1) = y(1);
    x_vel(1) = x_vel(find(~failpoints,1,'first'));
    y_vel(1) = y_vel(find(~failpoints,1,'first'));
    failpoints(1) = false;
  end
  if failpoints(end)
    x_smooth(end) = x(end);
    y_smooth(end) = y(end);
    x_vel(end) = x_vel(find(~failpoints,1,'last'));
    y_vel(end) = y_vel(find(~failpoints,1,'last'));
    failpoints(end) = false;
  end
  x_smooth(failpoints) = interp1(t(~failpoints),x_smooth(~failpoints), ...
      t(failpoints),'spline');
  y_smooth(failpoints) = interp1(t(~failpoints),y_smooth(~failpoints), ...
      t(failpoints),'spline');
  x_vel(failpoints) = interp1(t(~failpoints),x_vel(~failpoints), ...
      t(failpoints),'linear');
  y_vel(failpoints) = interp1(t(~failpoints),y_vel(~failpoints), ...
      t(failpoints),'linear');
  if any(isnan(x_smooth)) || any(isnan(x_vel)) || ...
      any(isnan(y_smooth)) || any(isnan(y_vel))
    error('there is a mistake in this code');
  end
end % end function LOESS_TRAJECTORY


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOESS ITERATIVE SMOOTHER TO OBTAIN SMOOTHED HEAD DIRECTION ESTIMATE
% Censor head direction when front and back lights are too close together
function theta = loess_head_direction(t,xfront,yfront,xback,yback, ...
    halfwidth,num_iterations,min_front_back_separation);

  if any(~isfinite(xfront)) || any(~isfinite(yfront)) || ...
      any(~isfinite(xback)) || any(~isfinite(yback)) || ...
      any(~isfinite(t)) || any(diff(t) <= 0)
    error(['input vectors must be finite non-NaN, and t vector must be ' ...
        'monotonically increasing']);
  end
  % Represent the separation between the front and back LEDs as a complex
  % number. Note that when the front and back LEDs are poorly separated (e.g.
  % when the head direction estimate is noisy), z is close to zero.
  z = complex(xfront-xback,yfront-yback);
  % Exclude samples in which the front and back LEDs do not achieve separation
  % to accurately resolve head direction
  selectionmask = (abs(z) >= min_front_back_separation);
  z_selected = z(selectionmask);
  t_selected = t(selectionmask);
  n = nnz(selectionmask);
  theta_selected = nan([n 1]);
  % Degree of the local smoothing polynominal; quadratic is adequate
  LOESS_DEGREE = 2;
  % Turn off warning in the LSCOV function, as otherwise MATLAB spits out a
  % flood of warning statements in the console
  warning off MATLAB:lscov:RankDefDesignMat
  % Index look-up for nearby samples that fall within the smoothing bandwidth
  % (these are indices into t_selected and z_selected)
  i_nearby = find_nearby_mex(t_selected,t_selected,-halfwidth,+halfwidth);
  % Keep track of points where the loess fit fails, usually because of missing
  % samples in the neighborhood (this is a bitmask for z_selected)
  failpoints = false([n 1]);
  % Perform iterative loess at each position sample
  z_coeffs = zeros([n LOESS_DEGREE+1]);
  for i = 1:n
    % Construct degree-d Vandermonde matrix: first column is dth-order
    % coefficient, second column is (d-1)th order coefficient, ... dth column is
    % linear term coefficient, (d+1)th column is constant term (this is the
    % convention for POLYVAL, POLYDER, POLYFIT, etc.)
    lags = t_selected(i_nearby{i}) - t_selected(i); 
    basis = ones([numel(lags) LOESS_DEGREE+1]);
    for d = LOESS_DEGREE:-1:1
      basis(:,d) = lags .* basis(:,d+1);
    end
    % Grab observations in this window
    z_obs = z_selected(i_nearby{i});
    % Initialize regression weights by tricube weighting over the smoothing
    % bandwidth
    weights = tricube(lags/halfwidth);
    iter_count = 0;
    while true
      % Fit polynomial to local observations (weighted by weights vector)
      try
        z_coeffs(i,:) = lscov(basis,z_obs,weights);
      catch
        failpoints(i) = true;
        break;
      end
      if (iter_count < num_iterations)
        % Adjust the weights acccording to the residuals: normalize residuals
        % w.r.t. 6*median(residuals), then update weights by multiplying with
        % bi-square function of normalized residuals
        residuals = abs(polyval(z_coeffs(i,:),lags) - z_obs);
        weights = weights .* ...
            bisquare(residuals/(6*median(abs(residuals))));
        iter_count = iter_count + 1;
      else
        break;
      end
    end
  end
  theta_selected(~failpoints) = angle(z_coeffs(find(~failpoints),end));
  % Patch with interpolated values at points where polynomial fit failed
  theta_selected(failpoints) = interp1_angle(t_selected(~failpoints), ...
      theta_selected(~failpoints),t_selected(failpoints),pi);
  % Censor samples in which the front and back LEDs do not achieve separation
  % to accurately resolve head direction
  theta_selected(abs(z_coeffs(:,end)) < min_front_back_separation) = NaN;
  % Copy selected values to the output vector
  theta = nan(size(t));
  theta(selectionmask) = theta_selected;
end % end function LOESS_HEAD_DIRECTION


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPEATED RUNNING MEDIAN DETECTION OF STOPS
function stopped = rrm_arrests(t,x,y,halfwidths,min_stop_duration,epsilon)
  n = numel(t);
  % Sequentially compute running medians over x and y, iterating through the
  % halfwidths specified in params.rrm_halfwidths
  x_median = nan(size(x));
  y_median = nan(size(y));
  for hw = halfwidths
    % i_nearby{i} contains indices to all samples that are within halfwidth of
    % the ith sample
    i_nearby = find_nearby_mex(t,t,-hw,+hw);
    % Surprisingly, this for loop is much faster than arrayfun
    for i = 1:n
      x_median(i) = median(x(i_nearby{i}));
      y_median(i) = median(y(i_nearby{i}));
    end
    % Replace x and y with their respective running medians
    x = x_median;
    y = y_median;
  end
  % After applying the running median smoothers, identify continuous intervals
  % that satisfy the criteria min_stop_duration (time) and epsilon (space). For
  % each sample, look ahead to min_stop_duration later and determine whether all
  % samples in this range fall within epsilon
  stopped = false([n 1]);
  for i = 1:(n-1)
    i_range = i : ...
        (i + find(t((i+1):end) - t(i) >= min_stop_duration,1,'first'));
    if ~isempty(i_range)
      stopped(i_range) = stopped(i_range) | ...
          all(pdist([x(i_range) y(i_range)]) < epsilon);
    end
  end
end % end function RRM_ARRESTS


