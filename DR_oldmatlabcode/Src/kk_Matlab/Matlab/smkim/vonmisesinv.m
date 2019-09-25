function x = vonmisesinv(p,kappa)
%VONMISESINV Inverse of the von Mises cumulative distribution function (cdf)
%
%   X = VONMISESINV(P,KAPPA) returns the inverse cdf for the von Mises
%   distribution with mean zero and concentration parameter KAPPA, evaluated at
%   values in P. Because the von Mises distribution is circular, the cdf is
%   defined on the arbitrary support [-pi,+pi], and values within X must fall
%   within this support interval. Thus, the following conditions are imposed:
%
%   VONMISESINV(0,KAPPA) == -pi
%   VONMISESINV(0.5,KAPPA) == 0
%   VONMISESINV(1,KAPPA) == +pi
%
%   The size of X is the common size of P and KAPPA. A scalar input is
%   interpreted as an array filled with a single value to match the size of the
%   other input. X contains NaN values where P is NaN.
%
%   Reference:
%
%   Hill G.W. (1977) Algorithm 518. Incomplete Bessel function I0: the von Mises
%   distribution [S14]. ACM Transactions on Mathematical Software 3:279-284.
%
%   section 3.3.6 in Fisher N.I. (1993) _Statistical analysis of circular data._
%   Cambridge University Press.
%
%   See also VONMISEPDF, VONMISESCDF, VONMISESRND, VONMISESFIT
%
%Depends on:
%   VONMISESCDF (written by SMK)
%   NORMINV (MATLAB Statistics Toolbox)
%
%Written by SMK 2009 December 07.
%

error(['This function exhibits numerical instability and is disabled ' ...
    'until bugs are fixed']);

%{

if (exist('vonmisescdf') ~= 2)
  error('VONMISESINV depends on m-file VONMISESCDF (written by SMK)');
end
if (exist('norminv') ~= 2)
  error('VONMISESINV depends on m-file NORMINV (written by SMK)');
end

if ~isfloat(p) || ~isreal(p) || any(p(:) > 1) || any(p(:) < 0)
  error('P must be real floating-point between zero and one');
end
if ~isfloat(kappa) || ~isreal(kappa) || ~all(isfinite(kappa(:))) || ...
    any(isnan(kappa(:))) || any(kappa(:) <= 0)
  error('KAPPA must be real positive finite non-NaN floating-point');
end

% Check for size consistency
try
  test = vonmisescdf(zeros(size(p)),kappa);
catch
  error('Non-scalar arguments must match in size.');
end
sz = size(test);

% Iterative approximation described by Fisher (1992). This converges very slowly
% for kappa > 50. We speed it up by seeding with a normal approximation.
x = norminv(p,zeros(sz),sqrt(1./kappa));
x(x <= -pi) = -pi + eps;
x(x >= +pi) = +pi - eps;
c = log(besseli(0,kappa));
TOL = 1e-6;
MAX_ITER = 1e4;
d = inf(sz);
count = 0;
while any(abs(d(:)) > TOL)
  g = vonmisescdf(x,kappa) - p;
  d = exp(log(abs(g)) + c - kappa .* cos(x));
  x = x - sign(g) .* d;
  % Correct for round-off errors
  x(x < -pi) = -pi;
  x(x > +pi) = +pi;
  count = count + 1;
  %disp([count, max(abs(d(:)))]);
  if (count > MAX_ITER)
    error(['Estimate did not converge to tolerance (%f) after maximum ' ...
        'number of iterations (%d). Kappa may be too large'], ...
        TOL,MAX_ITER);
  end
end

%}

