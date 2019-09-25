function y = vonmisespdf(x,mu,kappa)
%VONMISESPDF von Mises probability density function (pdf)
%   Y = VONMISESPDF(X,MU,KAPPA) returns the pdf of the 
%   von Mises distribution with mean MU and concentration 
%   parameter KAPPA, evaluated at the values in X. The
%   size of Y is the common size of the input arguments. A 
%   scalar input is interpreted as a matrix filled with a 
%   single value to match the size of the other inputs.
%
%   See also VONMISESCDF, VONMISESINV, VONMISESRND, VONMISESFIT.
%
%Written by smk, 3 November 2008.
%

if ~isfloat(x) || ~isreal(x) || any(isinf(x(:)))
  error('X must be real non-Inf floating-point');
end
if any((x(:) < -pi) | (x(:) > +pi))
  %warning('Some elements of X are outside the bounds [-pi,+pi]');
end
if ~isfloat(mu) || ~isreal(mu) || ~all(isfinite(mu(:))) || any(isnan(mu(:)))
  error('MU must be real finite non-NaN floating-point');
end
if any((mu(:) < -pi) | (mu(:) > +pi))
  error('Some elements of MU are outside the bounds [-pi,+pi]');
end
if ~isfloat(kappa) || ~isreal(kappa) || ~all(isfinite(kappa(:))) || ...
    any(isnan(kappa(:))) || any(kappa(:) < 0)
  error('KAPPA must be real positive finite non-NaN floating-point');
end

try
  y = exp(kappa .* cos(x - mu)) ./ (2 * pi * besseli(0,kappa));
catch
  error('Non-scalar arguments must match in size.');
end


