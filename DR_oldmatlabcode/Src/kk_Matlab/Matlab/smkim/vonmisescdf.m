function p = vonmisescdf(x,kappa)
%VONMISESCDF Von Mises cumulative distribution function (cdf).
%
%   P = VONMISESCDF(X,KAPPA) returns the cdf of the von Mises distribution with
%   mean zero and concentration parameter KAPPA, evaluated at values in X.
%   Because the von Mises distribution is circular, the cdf is defined on the
%   arbitrary support [-pi,+pi], and values within X must fall within this
%   support interval. Thus, the following conditions are imposed:
%
%   VONMISESCDF(-pi,KAPPA) == 0
%   VONMISESCDF(0,KAPPA) == 0.5
%   VONMISESCDf(+pi,KAPPA) == 1
%
%   The size of P is the common size of X and KAPPA. A scalar input is
%   interpreted as an array filled with a single value to match the size of the
%   other input. P contains NaN values where X is NaN.
%
%   Reference:
%
%   Hill G.W. (1977) Algorithm 518. Incomplete Bessel function I0: the von Mises
%   distribution [S14]. ACM Transactions on Mathematical Software 3:279-284.
%
%   See also VONMISEPDF, VONMISESINV, VONMISESRND, VONMISESFIT.
%
%Depends on:
%   VONMISESPDF (written by SMK)
%
%Written by SMK 2009 December 07.
%

if (exist('vonmisespdf') ~= 2)
  error('VONMISESCDF depends on m-file VONMISESPDF (written by SMK)');
end

if ~isfloat(x) || ~isreal(x) || any(x(:) < -pi) || any(x(:) > +pi)
  error('X must be real floating-point on [-pi,+pi]');
end
if ~isfloat(kappa) || ~isreal(kappa) || ~all(isfinite(kappa(:))) || ...
    any(isnan(kappa(:))) || any(kappa(:) <= 0)
  error('KAPPA must be real positive finite non-NaN floating-point');
end

% Check for size consistency
try
  test = vonmisespdf(x,0,kappa);
catch
  error('Non-scalar arguments must match in size.');
end
sz = size(test);
p = nan(sz);
if isscalar(x)
  x = repmat(x,sz);
end
if isscalar(kappa)
  kappa = repmat(kappa,sz);
end
assert(isequal(size(p),size(kappa),size(x)));

% The zero-mean von Mises cdf can be approximated by a nested series expansion
% containing modified Bessel function ratios. The number of nested terms is
% increased until the approximation does not change more than TOL. The number
% of terms required is larger for larger kappa.
TOL = 1e-9;
idx = find(isnan(p) & (abs(x) < pi));
p_last = 0.5*ones(size(x(idx)));
% Start with 50 nested terms, increase as necessary
n = 50;
while true
  V = 0;
  for i = n:-1:1
    V = (besseli(i,kappa(idx)) ./ besseli(i-1,kappa(idx))) .* ...
        (sin(i*x(idx))/i + V);
  end
  p(idx) = 0.5 + x(idx)/(2*pi) + V/pi;
  if any(abs(p(idx) - p_last) > TOL)
    p_last = p(idx);
    n = n + 10;
  else
    break;
  end
end

%TODO: for large kappa, a Gaussian approximation may have better numerical
%precision.

% Due to round-off error, P may be slightly greater than one or slightly less
% than zero. We replace these by the limit values 0 or 1.
p(p > 1) = 1;
p(p < 0) = 0;
% Manually assign the boundary conditions to correct any round-off error.
p(x == -pi) = 0;
p(x == +pi) = 1;


