% RANDCOVA:  Randomized data for an ANCOVA model, allowing for unequal slopes.
%
%     Syntax: [x,y,grps] = randcova(n,[xmin xmax],regr_params,{residstd},{dp})
%
%           n = vector of sample sizes, one per group.
%           [xmin xmax] = vector of limits of independent variable.
%           regr_params = matrix of regression parameters: rows = [b0 b1],
%                           one row per group.
%           residstd =    standard deviation of residuals [default = 1].
%           dp =          number of decimal positions for resultant variables
%                           [default = 2].
%           ------------------------------------------------------------------
%           x,y =         column vectors of independent & dependent variables.
%           grps =        column vector of group identifiers (1,2,...).
%

% RE Strauss, 11/5/96
%   9/19/99 - updated handling of default input arguments.
%  11/29/99 - changed sequence of output arguments.

function [x,y,grps] = randcova(n,xbnd,params,residstd,dp)
  if (nargin < 4) residstd = []; end;
  if (nargin < 5) dp = []; end;

  if (isempty(residstd))
    residstd = 1;
  end;
  if (isempty(dp))
    dp = 2;
  end;
  dp = 10^dp;

  ngrps = length(n);                      % Number of groups
  nobs =  sum(n);                         % Number of observations

  xmin = xbnd(1);                         % Uniform-random independent variable
  xmax = xbnd(2);
  x = rand(nobs,1)*(xmax-xmin) + xmin;

  grps = zeros(nobs,1);
  y = zeros(nobs,1);

  first = 1;
  for g = 1:ngrps
    last = first-1 + n(g);
    grps(first:last) = g*ones(n(g),1);     % Group identifiers
    e = randn(n(g),1) * residstd;         % Residuals
    y(first:last) = params(g,1) + params(g,2)*x(grps==g) + e;
    first = last+1;
  end;

  x = round(x*dp)/dp;                   % Specified number of decimal positions
  y = round(y*dp)/dp;

  xrng = xmax - xmin;
  xmin = xmin - xrng/20;
  xmax = xmax + xrng/20;

  clf;
  plot(x,y,'oy');
  hold on;
  for g = 1:ngrps
    y1 = params(g,1)+params(g,2)*xmin;
    y2 = params(g,1)+params(g,2)*xmax;
    plot([xmin xmax],[y1 y2],'r');
  end;
  hold off;

  return;
