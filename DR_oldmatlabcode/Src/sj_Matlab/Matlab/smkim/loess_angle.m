function [b, r] = loess_angle(x,a,xi,params)
%LOESS_ANGLE Robust iterative locally-weighted linear-circular regression
%
%   B = LOESS_ANGLE(X,A,XI,PARAMS) returns the matrix B of loess coefficients
%   that best fit the observations in X and A at the specified grid points XI.
%   X is an N-by-P design matrix, with rows corresponding to observations and
%   columns to predictor variables. A is an N-by-1 vector of circular response
%   observations, measured in radians on the interval [-pi,+pi]. Values in A
%   that lie outside of this interval are wrapped to the equivalent angle on
%   this interval. XI is a NI-by-P matrix, with rows corresponding to
%   P-dimensional points in the domain of X. X, A, and XI must all be real and
%   floating-point.
%
%   LOESS_ANGLE performs locally-weighted linear-circular regression at each
%   point in XI. B is a NI-by-(P+1) array of the resulting regression
%   coefficients. By default, LOESS_ANGLE adds a column of ones to X,
%   corresponding to a constant term in the first column of B. Do not enter a
%   column of ones directly into the X matrix.
%
%   PARAMS is a scalar struct with the following fields:
%     'Method': a linear-circular regression method argument accepted by
%         LINEAR_CIRCULAR_FI 'maximum_concentration' is recommended.
%     'Shrinkage': the ridge parameter, a non-negative real finite scalar. Zero
%         means that there is no shrinkage penalty. Larger values more
%         aggressively penalize overfitting, with the tradeoff of increasing
%         bias in the estimated coefficients.
%     'DistFcn': specifies the distance metric of the local regression
%         bandwidth. This can be a function handle or a string argument of the
%         type that is accepted by PDIST (MATLAB Statistics Toolbox).
%     'WidthType': a string, either 'fixed' or 'proportional'. 'fixed'
%         means that the smoothing bandwidth is expressed in units of X, while
%         'proportional' means that the smoothing bandwidth is expressed as a
%         fraction of the data points.
%     'Width': a scalar smoothing parameter. If 'WidthType' is 'fixed', then
%         this parameter must be expressed in the same units as
%         PDIST(X,PARAMS.DistFcn); for each point in XI, all data points that
%         are within PARAM.Width contribute to the local regression. If
%         'WidthType' is 'proportional', then this parameter must be a real
%         floating-point value greater than zero and less than or equal to one;
%         the nearest PARAMS.Width fraction of data points to each point in XI
%         are contribute to the local regression.
%     'NumIter': a non-negative real integer which specifies the number
%         of re-weighting iterations for robust smoothing. On each iteration,
%         observations with large residuals are down-weighted. A value of zero
%         means that there should be no iterative reweighting.
%
%   [B,R] = LOESS_ANGLE(X,A,XI,PARAMS) returns a N-by-1 vector of angular
%   residuals corresponding to the N observations. Residuals are bounded on the
%   interval [-pi,+pi] due to the circular nature of the response.
%
%   LOESS_ANGLE(X,A,[],PARAMS) is equivalent to LOESS_ANGLE(X,A,X,PARAMS).
%
%   LOESS_ANGLE treats NaNs in X or A as missing values, and removes them.
%
%   See also: LOESS (written by SMK), LINEAR_CIRCULAR_FIT (written by SMK)
%
%   References:
%   [1] Cleveland W.S. (1979) Robust locally weighted regression and smoothing
%       scatterplots. _Journal of the American Statistical Association_
%       74:829-836.
%   [2] Rowe D.B., Meller C.P., Hoffmann R.G. (2007) Characterizing phase-only
%       fMRI data with an angular regression model. _Journal of Neuroscience
%       Methods_ 161:331-341.
%   [3] Schmidt R., Diba K., Leibold C., Schmitz D., Buzsaki G., Kempter R.
%       (2009) Single-trial phase precession in the hippocampus. _Journal of
%       Neuroscience_ 29:13232-13241.
%
%Depends on:
%   LINEAR_CIRCULAR_FIT (written by SMK)
%   PDIST (MATLAB Statistics Toolbox)
%   TRICUBE (MATLAB Statistics Toolbox)
%   BISQUARE (MATLAB Statistics Toolbox)
%
%Written by SMK, 2009 December 14.
%

% Dependencies
if (exist('linear_circular_fit') ~= 2)
  error(['LOESS_ANGLE depends on the m-file LINEAR_CIRCULAR_FIT ' ...
      '(written by SMK)']);
end
if (exist('pdist') ~= 2)
  error(['LOESS_ANGLE depends on the m-file PDIST ' ...
      '(MATLAB Statistics Toolbox)']);
end
if (exist('tricube') ~= 2)
  error(['LOESS_ANGLE depends on the m-file TRICUBE ' ...
      '(written by SMK)']);
end
if (exist('bisquare') ~= 2)
  error(['LOESS_ANGLE depends on the m-file BISQUARE ' ...
      '(written by SMK)']);
end

% Input checking
if ~isfloat(x) || ~isreal(x) || (ndims(x) > 2) || any(isinf(x(:)))
  error('X must be a matrix of real non-Inf floating-point values');
end
if ~isfloat(a) || ~isreal(a) || ~isvector(a) || any(isinf(a(:))) || ...
    (numel(a) ~= size(x,1))
  error(['A must be a vector of real non-Inf floating-point values ' ...
      'whose number of elements equals the number of rows in X']);
end
if any((a < -pi) | (a > +pi))
  warning('Some elements of A lie outside the interval [-pi,+pi]');
end
if ~isequal(xi,[]) && (~isfloat(xi) || ~isreal(xi) || (ndims(xi) > 2) || ...
    (size(xi,2) ~= size(x,2)) || ~all(isfinite(xi(:))))
  error(['XI must be either the empty array [] or a matrix of real ' ...
      'non-Inf, non-NaN floating-point values, with the same number of ' ...
      'columns as X']);
end
REQUIRED_PARAMS_FIELDS = { ...
    'Method', ...
    'Shrinkage', ...
    'DistFcn', ...
    'WidthType', ...
    'Width', ...
    'NumIter' };
if ~isstruct(params) || ~isscalar(params) || ...
    ~all(isfield(params,REQUIRED_PARAMS_FIELDS))
  error('PARAMS is not a scalar struct or is missing required fields');
end
if ~isscalar(params.Shrinkage) || ~isreal(params.Shrinkage) || ...
    ~isfinite(params.Shrinkage) || (params.Shrinkage < 0)
  error('''Shrinkage'' field of PARAMS must be a non-negative real scalar');
end
switch (params.WidthType)
case 'fixed'
  if ~isfloat(params.Width) || ~isscalar(params.Width) || ...
      ~isreal(params.Width) || (params.Width <= 0)
    error('params.Width must be a positive floating-point real scalar');
  end
case 'proportional'
  if ~isfloat(params.Width) || ~isscalar(params.Width) || ...
      ~isreal(params.Width) || (params.Width <= 0) || (params.Width > 1)
    error(['params.Width must be a floating-point real scalar ' ...
        'greater than zero and less than or equal to one']);
  end
otherwise
  error(['''WidthType'' field of PARAMS must be either ''fixed'' or ' ...
      '''proportional''']);
end
if ~isscalar(params.NumIter) || ~isreal(params.NumIter) || ...
    ~isfinite(params.NumIter) || (params.NumIter < 0) || ...
    (params.NumIter ~= round(params.NumIter))
  error('''NumIter'' field of PARAMS must be a non-negative integer');
end

% Allocate outputs
ni = size(xi,1);
p = size(x,2);
if ~isempty(xi)
  b = nan([ni, 1+p]);
else
  b = nan([size(x,1) 1+p]);
end
r = nan(size(a));

% For convenience, reshape A to be a column vector
if (size(a,1) < size(a,2))
  a = a';
end

% Exclude NaNs. If there are no NaNs, then valid_obs == 1:numel(a).
valid_obs = find(~any(isnan(x),2) & ~any(isnan(a)));
n = numel(valid_obs);
% We are changing X and A. Don't be confused by discrepancies with the original
% inputs.
x = x(valid_obs,:);
a = a(valid_obs);

% Return NaN outputs if there are no valid data points
if isempty(a)
  return;
end

% Test that DistFcn works on the concatenated matrix [X; XI]
try
  d = squareform(pdist([x; xi],params.DistFcn));
  assert( isfloat(d) && isreal(d) && ~any(isnan(d(:))) && ~any(d(:) < 0) );
catch
  error(['''DistFcn'' field of PARAMS does not successfully return ' ...
      'real floating-point non-NaN non-negative distances']);
end
% Remove columns to produce an N-by-(N+NI) rectangular matrix. Each row
% corresponds to a row of the original input X, and each column corresponds to
% a row in the concatenated matrix [X; XI]
d = d(1:n,:);

% For each row in [X; XI], find nearby rows in X
i_local = cell([n+ni, 1]);
switch params.WidthType
case 'fixed'
  for i = 1:(n+ni)
    i_local{i} = find(abs(d(:,i)) < params.Width);
  end
case 'proportional'
  n_obs = ceil(params.Width * n);
  for i = 1:(n+ni)
    [vals, j] = sort(abs(d(:,i)),1,'ascend');
    i_local{i} = j(1:n_obs);
  end
otherwise
  error('There is an error in LOESS_ANGLE');
end

% Perform robust iterative linear-circular fit at each point in [X; XI]
for i = 1:(n+ni)
  if (i <= n)
    % We only need to fit at points in X if the caller requested residuals or
    % if no XI was specified
    if ~isempty(xi) && (nargout < 2)
      continue;
    end
    x_target = x(i,:);
  else % if (i > n)
    % This is a point in XI
    x_target = xi(i-n,:);
  end
  % Center the predictors around x_target
  predictors = bsxfun(@minus,x(i_local{i},:),x_target);
  responses = a(i_local{i});
  % Initial tri-cube weighting
  d_local = d(i_local{i},i);
  switch params.WidthType
  case 'fixed'
    weights = tricube(d_local ./ params.Width);
  case 'proportional'
    weights = tricube(d_local ./ max(d_local));
  end
  % We estimate a (1+P)-element vector of regression coefficients at each point
  coeffs = zeros([1, 1+p]);
  try
    coeffs = linear_circular_fit(predictors,responses,weights, ...
        params.Method,params.Shrinkage);
  catch
    warning(['Could not estimate fit within local bandwidth %f ', ...
        'of %s'],params.Width,mat2str(x_target));
    continue;
  end
  count = 0;
  while (params.NumIter > count)
    % Normalize residuals (bounded on the interval [-pi,+pi]) of the local
    % observations with respect to 7*std_angle(residuals), then reweight by a
    % bisquare function of the normalized residuals
    residuals = minus_angle(responses, ...
        coeffs(1) - (predictors * coeffs(2:end)'));
    weights = weights .* bisquare(residuals/(7*std_angle(residuals)));
    % Refit
    try
      coeffs = linear_circular_fit(predictors,responses,weights, ...
          params.Method,params.Shrinkage);
    catch
      warning(['Could not estimate fit within local bandwidth %f ', ...
          'of %s'],params.Width,mat2str(x_target));
      continue;
    end
    count = count + 1;
  end
  % Assign to appropriate output. Remember that outputs match the sizes of the
  % original inputs, which may not be the same as the censored X and A; we need
  % to look up indices in valid_obs
  if (i <= n)
    r(valid_obs(i)) = a(i) - coeffs(1);
    if isempty(xi)
      % If XI is not specified, then estimate fit at each data point
      b(valid_obs(i),:) = coeffs;
    end 
  else % if (i > n)
    b(i-n,:) = coeffs;
  end
end


