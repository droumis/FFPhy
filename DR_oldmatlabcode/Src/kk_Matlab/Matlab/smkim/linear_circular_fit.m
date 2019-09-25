function [b, r] = linear_circular_fit(x,a,w,method,penalty)
%LINEAR_CIRCULAR_FIT Wrapped linear regression with circular response variable.
%
%   B = LINEAR_CIRCULAR_FIT(X,A,W,METHOD) returns the vector B of regression
%   coefficients, estimated to fit the wrapped linear model
%
%     A_HAT = (B(1) + B(2:end)*X + pi, 2*pi) - pi
%
%   X is a N-by-P design matrix, with rows corresponding to observationsa nd
%   columns to predictor variables. A is an N-by-1 vector of circular response
%   observations, with values in radians on the interval [-pi,+pi]. Values in A
%   that lie outside of this interval will be wrapped to the equivalent angle
%   (warning will be given). W is an N-element vector of real non-negative
%   finite scalar weights. B is a 1-by-(P+1) vector of regression coefficients.
%   By default, LINEAR_CIRCULAR_FIT adds a column of ones to X, corresponding
%   to a constant term in the first element of B. Do not enter a column of ones
%   directly into the X matrix.
%
%   METHOD can be one of the following strings:
%     'maximum_concentration': maximize the circular concentration of the
%         angular residuals
%     'Fisher_Lee': maximum-likelihood fit assuming homoscedastic von Mises
%         distribution of residuals and arctangent link function 
%         X:[-Inf,+Inf] -> A:[-pi,+pi] (Fisher & Lee, 1992).
%
%   B = LINEAR_CIRCULAR_FIT(X,A,W,METHOD,PENALTY) imposes a shrinkage penalty
%   on the sum of the squares of the fit coefficients (after standardizing X),
%   a.k.a. ridge regression. This is useful because linear-circular regression
%   can overfit the data with arbitrarily large coefficients. PENALTY must be a
%   non-negative finite real scalar.
%
%   LINEAR_CIRCULAR_FIT(X,A,[],METHOD,PENALTY) is equivalent to
%   LINEAR_CIRCULAR_FIT(X,A,ones(size(A)),METHOD,PENALTY).
%
%   References:
%   [1] Fisher N.I. & Lee A.J. (1992) Regression models for an angular
%       response. _Biometrics_ 48:665-677.
%   [2] Rowe D.B., Meller C.P., Hoffmann R.G. (2007) Characterizing phase-only
%       fMRI data with an angular regression model. _Journal of Neuroscience
%       Methods_ 161:331-341.
%   [3] Schmidt R., Diba K., Leibold C., Schmitz D., Buzsaki G., Kempter R.
%       (2009) Single-trial phase precession in the hippocampus. _Journal of
%       Neuroscience_ 29:13232-13241.
%
%Depends on:
%   RECENTER_ANGLE (written by SMK)
%   MINUS_ANGLE (written by SMK)
%   VONMISESFIT (written by SMK)
%   FMINCON (MATLAB Optimization Toolbox)
%   FMINUNC (MATLAB Optimization Toolbox)
%
%Written by SMK, 2009 December 11.
%

  % Check dependencies
  if (exist('minus_angle') ~= 2)
    error(['LINEAR_CIRCULAR_FIT depends on m-file MINUS_ANGLE ' ...
        '(written by SMK)']);
  end
  if (exist('vonmisesfit') ~= 2)
    error(['LINEAR_CIRCULAR_FIT depends on m-file VONMISESFIT ' ...
        '(written by SMK)']);
  end
  if (exist('fmincon') ~= 2)
    error(['LINEAR_CIRCULAR_FIT depends on m-file FMINCON ' ...
        '(MATLAB Optimization Toolbox)']);
  end
  if (exist('fminunc') ~= 2)
    error(['LINEAR_CIRCULAR_FIT depends on m-file FMINUNC ' ...
        '(MATLAB Optimization Toolbox)']);
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
  if isequal(w,[]) 
    w = ones(size(a));
  elseif ~isfloat(w) || ~isreal(w) || ~isvector(w) || ...
      (numel(w) ~= numel(a)) || ~all(isfinite(w)) || ~all(w >= 0)
    error(['W must be either the empty array [] or a vector of real ' ...
        'finite non-negative floating-point values, with the same number ' ...
        'of elements as A']);
  end
  switch method
  case 'maximum_concentration'
    fit_routine = @maximum_concentration;
  case 'Fisher_Lee'
    fit_routine = @Fisher_Lee;
  otherwise
    error('METHOD must be one of the recognized string codes');
  end
  if (nargin < 5)
    penalty = 0
  elseif ~isreal(penalty) || ~isscalar(penalty) || ~ isfinite(penalty) || ...
      (penalty < 0)
    error('PENALTY must be a non-negative real finite scalar');
  end

  % Allocate outputs
  p = size(x,2);
  b = nan([1, p+1]);
  r = nan(size(a));

  % For convenience, reshape A and W to be a column vector
  if (size(a,1) < size(a,2))
    a = a';
  end
  if (size(w,1) < size(w,2))
    w = w';
  end

  % Exclude NaNs. If there are no NaNs, then valid_obs == 1:numel(y).
  valid_obs = find(~any(isnan(x),2) & ~any(isnan(a)));
  n = numel(valid_obs);
  % We are changing X, A, and W. Don't be confused by discrepancies with the
  % original inputs.
  x = x(valid_obs,:);
  a = a(valid_obs);
  w = w(valid_obs);

  % Return NaN outputs if there are no valid data points
  if isempty(a)
    return;
  end

  % Use the rank-revealing QR to check for linearly-dependent columns in X.
  [Q,R,perm] = qr(x,0);
  % Skip if there are not enough data points
  COND_MAX = 1e+6;
  if (size(R,1) < size(R,2)) || (condest(R) > COND_MAX)
    warning('Not enough observations are available for well-conditioned fit');
    return;
  end
  p_reduced = sum(abs(diag(R)) > max(n,p)*eps(R(1)));
  if (p_reduced < p)
    warning('Basis for linear regression is rank deficient');
    perm = perm(1:p_reduced);
  end

  % Perform linear-circular regression on reduced set of X columns, skipping
  % columns that were thrown out because of linear dependence. The first
  % element of B is reserved for the constant offset term
  b = zeros([1, 1+p]);
  %try
    [b(1), b(1+perm)] = fit_routine(x(:,perm),a,w,penalty);
  %catch
  %  error('Error in linear-circular regression routine');
  %end

  % Remember that output R matches the size of the original input, which may
  % not be the same as the censored X and Y; we need to look up indices in
  % valid_obs
  r(valid_obs) = minus_angle(a(valid_obs),b(1));

end % end main function LINEAR_CIRCULAR_FIT

function [offset, coeffs] = maximum_concentration(x,a,w,penalty)
  % The linear-circular model is
  %     a_hat = offset + coeffs*x + VM(0,kappa)
  % and the objective is to maximize the concentration of the angular
  % residuals subject to a shrinkage penalty.
  %
  % (1) Standardize the columns of X to have zero mean and unity variance. This
  %     is necessary so that the shrinkage penalty is fair.
  % (2) Estimate B to maximize the resultant length of the weighted sum of the
  %     angular residual phasors.
  % (3) Given the estimated B, find constant term B0 that maximizes the real
  %     component of the weighted sum of the angular residual phasors.
  % (4) Scale B0 and B according to the mean and standard deviation of each
  %     column of X
  p = size(x,2);
  % Normalize the columns of X to zero mean and unity variance.
  x_mean = mean(x);
  x_std = std(x,1,1);
  % If standard deviation is very small or NaN, replace with unity
  x_std(isnan(x_std) | (x_std < sqrt(eps(class(x_std))))) = 1;
  x = bsxfun(@rdivide,bsxfun(@minus,x,x_mean),x_std);
  % For pair of angles phi and rho, abs(exp(sqrt(-1)*phi)/exp(sqrt(-1)*rho))
  % is within floating-point tolerance of unity, so we can take a weighted sum
  % of these complex phasor ratios to get something that is proportional to
  % mean resultant length.
  objective1 = @(arg) -abs(sum(w .* exp(sqrt(-1)*(a - arg*x)))) + ...
        penalty*(arg*arg');
  options1 = optimset( ...
      'LargeScale'    , 'off', ...
      'FunValCheck'   , 'on', ...
      'Diagnostics'   , 'off', ...
      'Display'       , 'notify', ...
      'TypicalX'      , zeros([1, p]), ...
      'FinDiffType'   , 'central' );
  coeffs = fminunc(objective1,zeros([1, p]),options1);
  residuals = a - coeffs*x;
  objective2 = @(arg) -sum(w .* cos(residuals - arg));
  options2 = optimset( ...
      'FunValCheck'   , 'on', ...
      'Diagnostics'   , 'off', ...
      'Display'       , 'off', ...
      'TypicalX'      , 0 , ...
      'FinDiffType'   , 'central' );
  mu = vonmisesfit(a);
  offset = fmincon(objective2,mu,[],[],[],[],-pi,+pi,[],options2);
  % Scale to the mean and standard deviation of X
  coeffs = coeffs ./ x_std;
  offset = minus_angle(offset, coeffs * x_mean');
end


function [offset, coeffs] = Fisher_Lee(x,a,w,penalty)
  % The Fisher-Lee model is
  %     a_hat ~ VM(offset + 2*atan(coeffs*x),kappa)
  % and the objective is to maximize the weighted log likelihood subject to
  % shrinkage penalty.
  n = size(x,1);
  p = size(x,2);
  % Normalize the columns of X to zero mean and unity variance.
  x_mean = mean(x);
  x_std = std(x,1,1);
  % If standard deviation is very small or NaN, replace with unity
  x_std(isnan(x_std) | (x_std < sqrt(eps(class(x_std))))) = 1;
  x = bsxfun(@rdivide,bsxfun(@minus,x,x_mean),x_std);
  % Negative weighted log-likelihood with shrinkage penalty. arg(1) is the
  % concentration parameter of the von Mises distribution, which we don't use
  % here. arg(2) is the constant term. arg(3:end) are the slope coefficients.
  objective = @(arg) n*log(besseli(0,arg(1))) - arg(1) * sum( ...
      w .* cos(a - arg(2) - 2*atan(arg(3:end)*x))) + ...
      penalty*(arg(3:end)*arg(3:end)');
  options = optimset( ...
      'LargeScale'    , 'off', ...
      'FunValCheck'   , 'on', ...
      'Diagnostics'   , 'off', ...
      'Display'       , 'notify', ...
      'TypicalX'      , zeros([1, 2+p]), ...
      'FinDiffType'   , 'central' );
  [mu0, kappa0] = vonmisesfit(a);
  estimates = fminunc(objective,[kappa0, mu0, zeros([1 p])], ...
      options);
  offset = estimates(2);
  % Correct for the multiplicative factor of 2 in the arctangent link function
  coeffs = 2*estimates(3:end)
  % Scale to the mean and standard deviation of X
  coeffs = coeffs ./ x_std;
  offset = minus_angle(offset, coeffs * x_mean');
end


