function [muhat, kappahat, muci, kappaci] = vonmisesfit(x,alpha)
%VONMISESFIT Parameter estimates and confidence intervals for von Mises circular data.
%   [MUHAT,KAPPAHAT] = VONMISESFIT(X) returns estimates of the parameters of the
%   von Mises distribution given the data in X.  MUHAT is an estimate of the
%   mean, and KAPPAHAT is an estimate of the concentration. If X is an array, 
%   then parameters are estimated for each (hyper)column along the first
%   dimension, and 
%
%     size(MUHAT) = size(KAPPAHAT) = [1, size(X,2), size(X,3), ...]
%
%   [MUHAT,KAPPAHAT,MUCI,KAPPACI] = VONMISESFIT(X) returns 95% confidence
%   intervals (computed by parametric bootstrap) for the parameter estimates.
%
%   [MUHAT,KAPPAHAT,MUCI,KAPPACI] = VONMISESFIT(X,ALPHA) returns 100(1-ALPHA)
%   percent confidence intervals for the parameter estimates.
%
%   References:
%   
%   Best D.J., Fisher N.I. (1981) The biaso f the maximum likelihood estimators
%   of the von Mises-Fisher concentration parameters. _Communications in
%   Statistics - Simulation and Computation_. B10:493-502.
%
%   Section 4.5.5. (pages 88-93) in Fisher N.I. (1993) _Statistical Analysis of
%   Circular Data_. Cambridge University Press.
%
%   See also: VONMISESPDF, VONMISESCDF, VONMISESINV, VONMISESRND.
%
%Depends on:
%   MINUS_ANGLE (written by SMK)
%   FSOLVE (MATLAB Optimization Toolbox)
%   VONMISESRND (written by SMK)
%   QUANTILE (MATLAB Statistics Toolbox)
%
%Written by SMK, 2009 December 9.
%
%TODO: support input that contains NaN elements!

if (exist('minus_angle') ~= 2)
  error(['VONMISESFIT depends on the m-file MINUS_ANGLE ' ...
      '(written by SMK)']);
end
if (exist('fsolve') ~= 2)
  error(['VONMISESFIT depends on the m-file FSOLVE ' ...
      '(MATLAB Optimization Toolbox)']);
end
if (exist('vonmisesrnd') ~= 2)
  error(['VONMISESFIT depends on the m-file VONMISESRND ' ...
      '(written by SMK)']);
end
if (exist('quantile') ~= 2)
  error(['VONMISESFIT depends on the m-file QUANTILE ' ...
      '(MATLAB Statistics Toolbox)']);
end

if ~isfloat(x) || isempty(x) || ~isreal(x) || any(isinf(x(:))) || ...
    any(isnan(x(:)))
  error('X must be an array of real non-NaN finite floating-point values');
end
if any((x(:) < -pi) | (x(:) > +pi))
  warning(['Some elements of X lie outside the interval [-pi,+pi]; ' ...
      'these will be wrapped']);
end
if isvector(x) && (size(x,1) < size(x,2))
  x = x';
end
n = size(x,1);

if (nargin < 2)
  alpha = 0.05;
end
if ~isscalar(alpha) || ~isreal(alpha) || ~isfloat(alpha) || ...
    ~((alpha < 1) && (alpha > 0))
  error(['ALPHA must be a real floating-point scalar that is ' ...
      'greater than zero and less than one']);
end

% Maximum likelihood estimate of mean is simply the circular mean
muhat = mean_angle(double(x),1);

% This piecewise approximation of the maximum likelihood estimation of the
% concentration parameter is given in Fisher (1993)
R_bar = abs(mean(complex(cos(double(x)),sin(double(x))),1));
kappa_seed = ...
    (2*R_bar + R_bar.^3 + 5/6*R_bar.^5) .* (R_bar < 0.53) + ...
    (-0.4 + 1.39*R_bar + 0.43./(1 - R_bar)) .* ...
    ((R_bar >= 0.53) & (R_bar < 0.85)) + ...
    (1./(R_bar.^3 - 4*R_bar.^2 + 3*R_bar)) .* (R_bar >= 0.85);
% Starting from this approximation, compute the exact maximum likelihood
% estimate by numerical procedure
try
  kappahat = fsolve(@(k) besseli(1,k) ./ besseli(0,k) - R_bar,kappa_seed, ...
      optimset('Display','off'));
catch
  error('Unable to compute maximum likelihood estimate of concentration');
end
% Bias correction for small sample size (Best & Fisher, 1981)
if (n <= 15)
  %warning(['Maximum likelihood estimate of kappa will be corrected for ' ...
  %    'small sample-size bias']);
  kappahat = max(kappahat - 2/(n * kappahat),0) .* (kappahat < 2) + ...
      (n-1)^3 * kappahat / (n^3 + n) .* (kappahat >= 2);
end

% Compute confidence intervals by parametric bootstrap: we draw samples from
% VM(muhat,kappahat) using VONMISEESNRD and estimate mu and kapp from each
% sample using (recursively) VONMISESFIT.
if (nargout > 2)
  if (n < 10)
    warning(['Due to small sample size (n = %d), bootstrap confidence ' ...
        'intervals may not be reliable']);
  end
  % Number of bootstrap simulations
  NSIM = 1e3;

  if isvector(x)
    assert(isscalar(muhat) && isscalar(kappahat));
    mu_sim = nan([NSIM 1]);
    kappa_sim = nan([NSIM 1]);
    for i = 1:NSIM
      r = vonmisernd(muhat,kappahat,size(x));
      [mu_sim(i), kappa_sim(i)] = vonmisesfit(r);
    end
  else
    sz = size(x);
    mu_sim = nan([NSIM, sz(2:end)]);
    kappa_sim = nan([NSIM, sz(2:end)]);
    % Cell array of subscripts into mu_sim and kappa_sim
    subs = arrayfun(@(j) 1:j,sz,'UniformOutput',false);
    repsz = [n, ones([1 ndims(x)-1])];
    for i = 1:NSIM
      r = vonmisesrnd(repmat(muhat,repsz),repmat(kappahat,repsz),size(x));
      subs{1} = i;
      [mu_sim(subs{:}) kappa_sim(subs{:})] = vonmisesfit(r);
    end
  end
  % Take quantiles of the bootstrap parameter estimates
  p_lower = alpha/2;
  p_upper = 1 - alpha/2;
  muci = [quantile(mu_sim,p_lower,1); quantile(mu_sim,p_upper,1)];
  kappaci = [quantile(kappa_sim,p_lower,1); quantile(kappa_sim,p_upper,1)];

end


