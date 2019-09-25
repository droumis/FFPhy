function r = vonmisesrnd(mu,kappa,sz)
%VONMISESRND Random arrays from the von Mises distribution
%
%   R = VONMISESRND(MU,KAPPA) returns an array of random angles chosen from the
%   von Mises distribution with mean MU and concentration parameter KAPPA. The
%   size of R is the common size of MU and KAPPA if both are arrays. If eitehr
%   parameter is a scalar, the size of R is the size of the other parameter.
%
%   R = VONMISESRND(MU,KAPPA,[M,N,...]) returns an M-by-N-by-... array.
%
%   Reference:
%
%   Best D.J., Fisher N.I. (1979) Efficient simulation of the von Mises
%   distribution. _Journal of the Royal Statistical Society. Series C (Applied
%   Statistics)_. 28: 152-157.
%
%   See also: VONMISESPDF, VONMISESCDF, VONMISESINV.
%
%Depends on:
%   SUM_ANGLE (written by SMK)
%   VONMISESPDF (written by SMK)
%   
%Written by SMK, 2009 December 7.
%

if (exist('sum_angle') ~= 2)
  error('VONMISESRND depends on m-file SUM_ANGLE (written by SMK)');
end
if (exist('vonmisespdf') ~= 2)
  error('VONMISESRND depends on m-file VONMISESPDF (written by SMK)');
end

if ~isfloat(mu) || ~isreal(mu) || ~all(isfinite(mu(:))) || any(isnan(mu(:)))
  error('MU must be real finite non-NaN floating-point');
end
if any((mu(:) < -pi) | (mu(:) > +pi))
  warning('Some elements of MU are outside the bounds [-pi,+pi]');
end
if ~isfloat(kappa) || ~isreal(kappa) || ~all(isfinite(kappa(:))) || ...
    any(isnan(kappa(:))) || any(kappa(:) < 0)
  error('KAPPA must be real non-negative finite non-NaN floating-point');
end

% Check for size consistency
if (nargin < 3)
  sz = [1 1];
end
try
  test = vonmisespdf(sum_angle(zeros(sz),mu),mu,kappa);
catch
  error('Non-scalar arguments must match in size.');
end
sz = size(test);
n = prod(sz);

% The algorithm by Best & Fisher (1979) uses a wrapped Cauchy density as an
% upper envelope for Monte Carlo importance sampling. The parameter of the
% wrapped Cauchy distribution is tuned to maximize the acceptance ratio.
if isscalar(kappa)
  kappa = repmat(kappa,sz);
end
tau = 1 + sqrt(1 + 4*kappa.^2);
rho = (tau - sqrt(2*tau)) ./ (2*kappa);

% Now generate random numbers from the wrapped Cauchy distribution and accept
% them at the correct proportion to obtain random numbers from a von Mises
% distribution with zero mean.
r = (1 + rho.^2) ./ (2*rho);
for i = 1:n
  % Draw from uniform on [-pi,+pi] if kappa equals zero
  if (kappa(i) == 0)
    r(i) = 2*pi*(rand(1) - 0.5);
    continue;
  end
  while true
    % We draw three random variables from uniform distribution on [0,1]
    u = rand([3 1]);
    % We use u(1) to seed the random draw from the wrapped Cauchy distribution
    z = cos(pi*u(1));
    f = (1 + z*r(i))/(z + r(i));
    % Ensure that f does not exceed unity (will lead to complex values)
    if (f > 1)
      f = 1;
    end
    c = kappa(i)*(r(i) - f);
    % We use u(2) to randomly accept/reject the desired proportion below the
    % upper envelope
    if (c*(2 - c) - u(2) > 0) || (log(c/u(2)) + 1 - c >= 0)
      % ACOS returns values on the interval [0, pi]; we use u(3) to randomly
      % flip the sign of half of the values
      r(i) = sign(u(3) - 0.5) * acos(f);
      break;
    end
  end
end

% Shift the values to be centered at MU and wrap within interval [-pi,+pi]
r = sum_angle(r,mu);


