function [a_smooth, slope] = linear_circular_loess(x,a,x_grid,params)
%LINEAR_CIRCULAR_LOESS Loess smoothing of circular response on linear predictor.
%
%   A_SMOOTH = LINEAR_CIRCULAR_LOESS(X,A,X_GRID,PARAMS) performs robust
%   iterative loess regression of a circular scalar response against a linear
%   scalar predictor. X is an N-element vector of real floating-point non-NaN
%   finite elements, A is an N-element vector of real floating-point non-NaN
%   finite elements (angles in radians), and X_GRID is a vector of neighborhood
%   centers at which the loess fit is to be estimated. PARAMS is a struct with
%   the following fields:
%     degree: degree of the local polynomial
%     bandwidth_type: either 'fixed' (an interval in the same units as X) or
%         'proportional' (as a fraction of all data points)
%     bandwidth: a scalar, either a number in the same units as X (if
%         bandwidth_type is 'fixed') or a fraction between 0 and 1 (if
%         bandwidth_type is 'proportional')
%     num_iterations: number of iterations to perform for robust fit, integer
%         greater than or equal to 1
%
%   [A_SMOOTH, SLOPE] = LINEAR_CIRCULAR_LOESS(...) also returns a vector of fit
%   slopes at the points in X_GRID, expressed in radians/[unit of X]
%
%   The regression model is a generalized form of the one proposed by Fisher &
%   Lee (1992), in which the circular values are mapped to the real continuum
%   by the tan function.
%
%   References:
%
%   Cleveland W.S. (1979) Robust locally weighted regression and smoothing
%   scatterplots. _Journal of the American Statistical Association_ 74:829-836.
%
%   Fisher N.I. & Lee A.J. (1992) Regression models for an angular response.
%   _Biometrics_ 48:665-677.
%
%Depends on:
%   FIND_NEARBY_MEX (written by SMK)
%   BISQUARE (written by SMK)
%   TRICUBE (written by SMK)
%   RECENTER_ANGLE (written by SMK)
%   MEAN_ANGLE (written by SMK)
%
%Written by SMK, 2009 November 20.
%

if (exist('find_nearby_mex') ~= 3)
  error(['LINEAR_CIRCULAR_LOESS depends on mex-file FIND_NEARBY_MEX ' ...
      '(written by smk)']);
end
if (exist('bisquare') ~= 2)
  error(['LINEAR_CIRCULAR_LOESS depends on m-file BISQUARE ' ...
      '(written by smk)']);
end
if (exist('tricube') ~= 2)
  error(['LINEAR_CIRCULAR_LOESS depends on m-file TRICUBE ' ...
      '(written by smk)']);
end
if (exist('recenter_angle') ~= 2)
  error(['LINEAR_CIRCULAR_LOESS depends on m-file RECENTER_ANGLE' ...
      '(written by smk)']);
end
if (exist('mean_angle') ~= 2)
  error(['LINEAR_CIRCULAR_LOESS depends on m-file MEAN_ANGLE' ...
      '(written by smk)']);
end

if ~isfloat(x) || ~isvector(x) || ~isreal(x) || ~all(isfinite(x))
  error('X must be a vector of real floating-point finite values');
end
if ~isfloat(x_grid) || ~isa(x_grid,class(x)) || ~isvector(x_grid) || ...
    ~isreal(x_grid) || ~all(isfinite(x_grid))
  error(['X_GRID must be a vector of real floating-point finite values ' ...
      'of the same numeric class as X']);
end
if ~isfloat(a) || ~isvector(a) || ~isreal(a) || ~all(isfinite(a))
  error('A must be a vector of real floating-point finite values');
end
if any((a > pi) | (a < pi))
  warning('some A values are outside the interval [-pi,+pi]');
end

REQUIRED_PARAMS = {
if ~isstruct(params) || ~isscalar(params) || ...
    ~all(isfield(params, ...
    {'degree','bandwidth_type','bandwidth','num_iterations'}))
  error('PARAMS is not valid');
end
if ~isscalar(params.degree) || ~isreal(params.degree) || ...
    ~(params.degree > 0) || ~(params.degree == round(params.degree))
  error('degree field of PARAMS must be a real positive integer');
end
if ~ischar(params.bandwidth_type) || ...
    ~any(strcmp(params.bandwidth_type,{'fixed','proportional'}))
  error(['bandwidth_type field of PARAMS must be either ' ...
      '''fixed'' or ''proportional''']);
end
switch params.bandwidth_type
case 'fixed'
  if ~isnumeric(params.bandwidth) || ~isscalar(params.bandwidth) || ...
      ~isreal(params.bandwidth) || ~isfinite(params.bandwidth) || ...
      ~(params.bandwidth > 0)
    error(['If bandwidth_type is ''fixed'', then bandwidth must be ' ...
        'a positive non-Nan real number']);
  end
case 'proportional'
  if ~isnumeric(params.bandwidth) || ~isscalar(params.bandwidth) || ...
      ~isreal(params.bandwidth) || ~isfinite(params.bandwidth) || ...
      ~(params.bandwidth > 0) || ~(params.bandwidth <= 1)
    error(['If bandwidth_type is ''proportional'', then bandwidth must be ' ...
        'a real scalar greater than zero and less than or equal to one']);
  end
otherwise
  error('%s is not a recognized value for bandwidth_type', ...
      params.bandwidth_type);
end
if ~isscalar(params.num_iterations) || ~isreal(params.num_iterations) || ...
    ~(params.num_iterations > 0) || ...
    ~(params.num_iterations == round(params.num_iterations))
  error('degree field of PARAMS must be a real positive integer');
end

% Depending on params.bandwidth_type, find nearest data points to each point in
% X_GRID
switch params.bandwidth_type
case 'fixed'
  i_nearby = find_nearby_mex(x_grid,x,-params.bandwidth/2,+params.bandwidth/2);
case 'proportional'
  % How many elements fit within the fractional bandwidth?
  n_nearby = ceil(numel(x)*params.bandwidth);
  i_nearby = cell([numel(x_grid) 1]);
  for i = 1:numel(x_grid)
    [d, j] = sort(abs(x_grid(i) - x),1,'ascend');
    i_nearby{i}] = j(1:n_nearby);
  end
otherwise
  error('something went wrong');
end

a_smooth = nan(size(x));
slope = nan(size(x));

% Turn off warning in the LSCOV function, as otherwise MATLAB spits out a
% flood of warning statements in the console
warning off MATLAB:lscov:RankDefDesignMat

% Now perform loess fit at each point in X_GRID
for i = 1:numel(x_grid)

  % Construct degree-d Vandermonde matrix: first column is dth-order
  % coefficient, second column is (d-1)th order coefficient, ... dth column is
  % linear term coefficient, (d+1)th column is constant term (this is the
  % convention for POLYVAL, POLYDER, POLYFIT, etc.)
  offsets = x(i_nearby{i}) - x_grid(i); 
  basis = ones([numel(offsets) params.degree+1]);
  for d = params.degree:-1:1
    basis(:,d) = offsets .* basis(:,d+1);
  end
  coeffs = zeros([n params.degree+1]);
  % Keep track of points where the fit fails. The fit can fail if X is poorly
  % spaced relative to X_GRID, or when params.bandwidth_type is 'fixed'
  failpoints = false([n 1]);
  % Initialize regression weights by tricube weighting over the smoothing
  % bandwidth
  switch params.bandwidth_type
  case 'fixed'
    weights = tricube(offsets/(params.bandwidth/2));
  case 'proportional'
    weights = tricube(offsets/max(abs(offsets)));
  end
  % Transform circular responses
  y = tan(recenter_angle(a,0)/2);

  iter_count = 0;
  while true
    % Fit polynomial to local observations (weighted by weights vector)
    try
      coeffs(i,:) = lscov(basis,x_obs,weights);
    catch
      % LSCOV will fail if there are not enough neighboring data points within
      % the bandwidth, in which case we flag this sample and move on to the
      % next sample
      failpoints(i) = true;
      break;
    end
    if (iter_count < params.num_iterations)
      % Adjust the weights acccording to the residuals: normalize residuals
      % w.r.t. 6*median(residuals), then update weights by multiplying with
      % bi-square function of normalized residuals
      
      %TODO
      %{
      residuals = polyval(coeffs(i,:),offsets) - x_obs;
      weights = weights .* ...
          bisquare(residuals/(6*median(abs(residuals))));
      iter_count = iter_count + 1;
      %}
    else
      break;
    end
  end
end



