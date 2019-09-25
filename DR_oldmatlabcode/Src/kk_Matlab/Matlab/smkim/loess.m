function [b, r] = loess(x,y,xi,params)
%LOESS Robust iterative locally-weighted linear regression
%
%   B = LOESS(X,Y,XI,PARAMS) returns the matrix B of loess coefficients that
%   best fit the observations in X and Y at the specified grid points XI. X is
%   an N-by-P design matrix, with rows corresponding to observations and columns
%   to predictor variables. Y is an N-by-1 vector of response observations. XI
%   is a NI-by-P matrix, with rows corresponding to P-dimensional points in the
%   domain of X. X, Y, and XI must all be real and floating-point.
%
%   LOESS performs locally-weighted linear regression at each point in XI. B is
%   a NI-by-(P+1) array of the resulting regression coefficients. By default,
%   LOESS adds a column of ones to X, corresponding to a constant term in the
%   first column of B. Do not enter a column of ones directly into the X matrix.
%   LOESS can also be used to perform locally-weighted polynomial regression by
%   constructing the appropriate matrix whose columns are polynomial expansions
%   of the predictors.
%
%   PARAMS is a scalar struct with the following fields:
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
%   [B,R] = LOESS(X,Y,XI,PARAMS) returns a N-by-1 vector of residuals
%   corresponding to the N observations.
%
%   LOESS(X,Y,[],PARAMS) is equivalent to LOESS(X,Y,X,PARAMS).
%
%   LOESS treats NaNs in X or Y as missing values, and removes them.
%
%   References:
%   [1] Cleveland W.S. (1979) Robust locally weighted regression and smoothing
%       scatterplots. _Journal of the American Statistical Association_
%       74:829-836.
%   [2] Cleveland W.S., Devlin S.J. (1988) Locally weighted regression: an
%       approach to regression analysis by local fitting. _Journal of the
%       American Statistical Association_ 83:596-610.
%
%Depends on:
%   PDIST (MATLAB Statistics Toolbox)
%   TRICUBE (written by SMK)
%   BISQUARE (written by SMK)
%
%Written by SMK, 2009 December 11.
%

% Dependencies
if (exist('pdist') ~= 2)
  error(['LOESS depends on the m-file PDIST ' ...
      '(MATLAB Statistics Toolbox)']);
end
if (exist('tricube') ~= 2)
  error(['LOESS depends on the m-file TRICUBE ' ...
      '(written by SMK)']);
end
if (exist('bisquare') ~= 2)
  error(['LOESS depends on the m-file BISQUARE ' ...
      '(written by SMK)']);
end

% Input checking
if ~isfloat(x) || ~isreal(x) || (ndims(x) > 2) || any(isinf(x(:)))
  error('X must be a matrix of real non-Inf floating-point values');
end
if ~isfloat(y) || ~isreal(y) || ~isvector(y) || any(isinf(y(:))) || ...
    (numel(y) ~= size(x,1))
  error(['Y must be a vector of real non-Inf floating-point values ' ...
      'whose number of elements equals the number of rows in X']);
end
if ~isequal(xi,[]) && (~isfloat(xi) || ~isreal(xi) || (ndims(xi) > 2) || ...
    (size(xi,2) ~= size(x,2)) || ~all(isfinite(xi(:))))
  error(['XI must be either the empty array [] or a matrix of real ' ...
      'non-Inf, non-NaN floating-point values, with the same number of ' ...
      'columns as X']);
end
REQUIRED_PARAMS_FIELDS = { ...
    'DistFcn', ...
    'WidthType', ...
    'Width', ...
    'NumIter' };
if ~isstruct(params) || ~isscalar(params) || ...
    ~all(isfield(params,REQUIRED_PARAMS_FIELDS))
  error('PARAMS is not a scalar struct or is missing required fields');
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
r = nan(size(y));

% For convenience, reshape Y to be a column vector
if (size(y,1) < size(y,2))
  y = y';
end

% Exclude NaNs. If there are no NaNs, then valid_obs == 1:numel(y).
valid_obs = find(~any(isnan(x),2) & ~any(isnan(y)));
n = numel(valid_obs);
% We are changing X and Y. Don't be confused by discrepancies with the original
% inputs.
x = x(valid_obs,:);
y = y(valid_obs);

% Return NaN outputs if there are no valid data points
if isempty(y)
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
  error('There is an error in LOESS');
end

% Criterion for determining whether a matrix is ill-conditioned
COND_MAX = 1e+6;

% Perform robust iterative fit at each point in [X; XI]
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
  responses = y(i_local{i});
  n_obs = size(predictors,1);
  % Append a leading column of ones for the constant term
  basis = [ ones([n_obs 1]), predictors ];
  % Use the rank-revealing QR to check for linearly-dependent columns
  [Q,R,perm] = qr(basis,0);
  % Skip if there are not enough data points within bandwidth; this can happen
  % if WidthType is 'fixed'
  if (size(R,1) < size(R,2)) || (condest(R) > COND_MAX)
    warning(['Not enough points are found within local bandwidth %f ', ...
        'of %s for well-conditioned fit'],params.Width,mat2str(x_target));
    continue;
  end
  p_reduced = sum(abs(diag(R)) > max(n_obs,p)*eps(R(1)));
  if (p_reduced < p)
    warning('Local basis for linear regression is rank deficient');
    perm = perm(1:p_reduced);
    % Always keep the first column of ones!
    if ~ismember(1,perm)
      perm = [1, perm(1:end-1)];
    end
  end
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
  % Perform weighted least-squares linear regression on reduced basis, skipping
  % basis components that were thrown out because of linear dependence
  try
    coeffs(perm) = linsolve(bsxfun(@times,basis(:,perm),weights), ...
        responses .* weights,struct('RECT',true));
  catch
    warning(['Could not estimate fit within local bandwidth %f ', ...
        'of %s'],params.Width,mat2str(x_target));
  end
  count = 0;
  while (params.NumIter > count)
    % Normalize residuals of the local observations with respect to
    % 6*median(abs(residuals)), then reweight by a bi-square function of the
    % normalized residuals
    residuals = responses - (basis * coeffs');
    weights = weights .* bisquare(residuals/(6*median(abs(residuals))));
    % Refit
    try
      coeffs(perm) = linsolve(bsxfun(@times,basis(:,perm),weights), ...
          responses .* weights,struct('RECT',true));
    catch
      warning(['Could not estimate fit within local bandwidth %f ', ...
          'of %s'],params.Width,mat2str(x_target));
    end
    count = count + 1;
  end
  % Assign to appropriate output. Remember that outputs match the sizes of the
  % original inputs, which may not be the same as the censored X and Y; we need
  % to look up indices in valid_obs
  if (i <= n)
    r(valid_obs(i)) = y(i) - coeffs(1);
    if isempty(xi)
      % If XI is not specified, then estimate fit at each data point
      b(valid_obs(i),:) = coeffs;
    end 
  else % if (i > n)
    b(valid_obs(i-n),:) = coeffs;
  end
end


