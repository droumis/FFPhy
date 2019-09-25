% Randmink: Random sample from a Minkowski distribution with parameter k>=0 
%           (k >= 0.2, in practice) and unit variance.  The distribution for 
%           k=2 is the standard normal distribution.
%
%     Usage: m = rankmink(k,r,c)  OR  m = randmink(k,{[r,c]})
%
%           k =   Minkowski parameter [recommend k > 0.2].
%           r,c = rows and columns of output matrix [default: r=1, c=1].
%           ------------------------------------------------------------
%           m =   [r x c] matrix of random values.
%

% RE Strauss, 12/25/99

function m = randmink(k,r,c)
  if (nargin < 1) k = []; end;
  if (nargin < 2) r = []; end;
  if (nargin < 3) c = []; end;

  if (isempty(k))
    error('  RANKMINK: parameter k required');
  else
    if (k < 0.25)
      disp('  RANKMINK warning: function impractical for k < 0.25');
    end;
  end;

  if (isempty(r))
    r = 1;
  end;
  if (isempty(c))
    if (length(r)>1)
      c = r(2);
      r = r(1);
    else
      c = 1;
    end;
  end;
  n = r*c;

  tol = 1e-6;                             % Tolerance for cdf bounds
  xn = 2500;                              % Number of points in empirical cdf

  xmax = 0;                               % Bounds for cdf
  fx = 1;
  while (fx > tol)
    xmax = xmax+1;
    fx = exp(-abs(xmax).^k);
  end;
  xmin = -xmax;

  x = linspace(xmin,xmax,xn);             % Empirical pdf
  fx = exp(-abs(x).^k);
  fx = fx./sum(fx);
  s = sqrt(sum(fx.*(x.^2)));              % Standard deviation of pdf
  fx = cumsum(fx);                        % Cdf

  u = rand(n,1);                          % Uniform-random values
  m = zeros(n,1);                         % Allocate output vector
  for i = 1:n                             % Find x indices corresponding to u
    m(i) = max(find(fx<=u(i)));
  end;
  m = (m./xn).*(xmax-xmin) + xmin;        % Transform to x values
  m = m ./ s;                             % Adjust to unit standard deviation

  if (c>1)                                % Reshape to form output matrix
    if (r==1)
      m = m';
    else
      m = reshape(m,r,c);
    end;
  end;

  return;

