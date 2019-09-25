function position = estimate_position(rawpos,varargin)
%ESTIMATE_POSITION Estimate head direction, smoothed position and velocity from raw video pixel coordinates of head-tracking LEDs.
%
%   POSITION = ESTIMATE_POSITION(RAWPOS,SMOOTHING_PARAMS) estimates
%   head direction, smoothed position and velocity from the coordinates of the
%   head-tracking LEDs in the RAWPOS.
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
%   Optional arguments:
%
%   centimeters_per_pixel
%     Parameter for scaling pixels. Default is 1 cm / pixel.
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
%     Default value is 2 cm.
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
%   BISQUARE (written by smk)
%   TRICUBE (writtten by smk)
%   INTERP1_ANGLE (writtten by smk)
%
%Written by smk (2009 November 1)
%
centimeters_per_pixel = 1;
front_back_marker_weights = [1 0];
loess_halfwidth = 0.75;
num_loess_iterations = 4;
rrm_halfwidths = [0.1 0.1 0.1];
epsilon = 0.5;
min_stop_duration = 0.5;
min_front_back_separation = 1;
TS_PER_SEC = 1e4;

otherArgs = procOptions(varargin);

if (exist('find_nearby_mex') ~= 3)
  error(['ESTIMATE_POSITION depends on mex-file FIND_NEARBY_MEX ' ...
      '(written by smk)']);
end
if (exist('bisquare') ~= 2)
  error('ESTIMATE_POSITION depends on m-file BISQUARE (written by smk)');
end
if (exist('tricube') ~= 2)
  error('ESTIMATE_POSITION depends on m-file TRICUBE (written by smk)');
end
if (exist('interp1_angle') ~= 2)
  error('ESTIMATE_POSITION depends on m-file INTERP1_ANGLE (written by smk)');
end

% Initialize POSITION and inherit meta-data from RAWPOS
position = cell(1,numel(rawpos));
params = struct('centimeters_per_pixel', centimeters_per_pixel, ...
  'front_back_marker_weights', front_back_marker_weights, ...
  'loess_halfwidth', loess_halfwidth, ...
  'num_loess_iterations', num_loess_iterations, ...
  'rrm_halfwidths', rrm_halfwidths, ...
  'epsilon', epsilon, ...
  'min_stop_duration', min_stop_duration, ...
  'min_front_back_separation', min_front_back_separation);
posinfo = struct('params',params, 'units', 'centimeters', ...
  'fields', 'time x y dir vel xvel yvel xres yres stopped');
[position{:}] = deal(posinfo);

for i = 1:numel(rawpos)
  if isempty(rawpos{i}.data)
    continue;
  end
  % Determine whether a back marker was tracked
  if (sum(sum(~isnan(rawpos{i}.data(:,4:5)))) == 0) || (sum(sum(rawpos{i}.data(:,4:5))) == 0)
    back_exists = false;
  else
    back_exists = true;
  end  
  % Convert timestamps to floating-point seconds
  t = double(rawpos{i}.data(:,1));
  % Promote to double to take advantage of MATLAB acceleration
  xfront = double(rawpos{i}.data(:,2));
  yfront = double(rawpos{i}.data(:,3));
  xback = double(rawpos{i}.data(:,4));
  yback = double(rawpos{i}.data(:,5));
  % Censor samples with NaN or zero values
  censor_mask = ...
      isnan(xfront) | isnan(yfront) | ((xfront == 0) & (yfront == 0));
      % isnan(xfront) | (xfront == 0) | isnan(yfront) | (yfront == 0);
  if (back_exists)
    censor_mask = censor_mask | ...
      isnan(xback) | isnan(yback) | ((xback == 0) & (yback == 0));
      % isnan(xback) | (xback == 0) | isnan(yback) | (yback == 0);
  end
  tuncensored = t;
  t(censor_mask) = [];
  xfront(censor_mask) = [];
  yfront(censor_mask) = [];
  if (back_exists)
    xback(censor_mask) = [];
    yback(censor_mask) = [];
  end
  % n is the number of elements in t *after* NaN samples are removed
  n = numel(t);
  xfront = xfront * centimeters_per_pixel;
  yfront = yfront * centimeters_per_pixel;
  if (back_exists)
    xback = xback * centimeters_per_pixel;
    yback = yback * centimeters_per_pixel;
  end

  % Inherit timestamps
  position{i}.data(:,1) = rawpos{i}.data(:,1);
  if (back_exists)
    % Estimate head position from weighted average of front and back
    if (size(front_back_marker_weights,1) > 1)
      % Coerce to be 2-element row vector for proper alignment
      front_back_marker_weights = front_back_marker_weights';
    end
    x = sum(bsxfun(@times,[xfront xback], ...
        front_back_marker_weights),2) ./ ...
        sum(front_back_marker_weights);
    y = sum(bsxfun(@times,[yfront yback], ...
        front_back_marker_weights),2) ./ ...
        sum(front_back_marker_weights);
  else
    x = xfront;
    y = yfront;
  end

  % Perform LOESS smoothing to get denoised position and velocity
  [x_smooth, y_smooth, x_vel, y_vel] = loess_trajectory( ...
      tuncensored,t,x,y, ...
      loess_halfwidth,num_loess_iterations);

  position{i}.data(:,2:3) = [x_smooth y_smooth];
  position{i}.data(:,5) = sqrt(sum([x_vel y_vel].^2,2)); % speed
  position{i}.data(:,6:7) = [x_vel y_vel]; % velocity
  position{i}.data(:,8:9) = nan;
  position{i}.data(~censor_mask,8) = x - x_smooth(~censor_mask);
  position{i}.data(~censor_mask,9) = y - y_smooth(~censor_mask);
  % Perform LOESS smoothing to get denoised head direction
  position{i}.data(:,4) = nan; % default head direction
  if (back_exists)
    position{i}.data(~censor_mask,4) = loess_head_direction( ...
      t,xfront,yfront,xback,yback, ...
      loess_halfwidth,num_loess_iterations, ...
      min_front_back_separation);
  end  
  % Find arrests by applying repeated running median filter
  %  to the smoothed trajectory
  position{i}.data(:,10) = rrm_arrests( ...
    tuncensored,x_smooth,y_smooth, ...
    rrm_halfwidths,min_stop_duration,epsilon);

  % disp(position(i));
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
  n = numel(tuncensored);
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
    lags = t(i_nearby{i}) - tuncensored(i); 
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
  x_smooth(failpoints) = interp1(tuncensored(~failpoints),x_smooth(~failpoints), ...
      tuncensored(failpoints),'spline');
  y_smooth(failpoints) = interp1(tuncensored(~failpoints),y_smooth(~failpoints), ...
      tuncensored(failpoints),'spline');
  x_vel(failpoints) = interp1(tuncensored(~failpoints),x_vel(~failpoints), ...
      tuncensored(failpoints),'linear');
  y_vel(failpoints) = interp1(tuncensored(~failpoints),y_vel(~failpoints), ...
      tuncensored(failpoints),'linear');
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
  if isempty(t_selected)
    theta = nan(length(t),1);
    return
  end
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
  % halfwidths specified in rrm_halfwidths
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


