% Testgrps: Constructs one or more multivariate-normal samples of 
%           equal variance and covariance (within sampling variation).
%           Centroids are selected randomly from uniform distributions. 
%
%     Usage: [x,g,m] = testgrps(k,n,p,{exact},{v},{r})
%
%           k =     number of groups to be created.
%           n =     sample size per group (if identical for all groups), 
%                     or vector (length k) of sample sizes.
%           p =     number of variables.
%           exact = optional boolean flag indicating, if true, that covariances 
%                     among random observations are to be exactly those of the 
%                     target matrix.  Default (=0): covariances vary randomly.
%           v =     optional within-group variance (identical for all 
%                     variables) [default = 2].
%           r =     optional within-group correlation (identical for all 
%                     pairs of variables) [default = 0.7].
%           ---------------------------------------------------------
%           x =     [k*n x p] data matrix.
%           g =     column vector (length k*n) of group identifiers 
%                     (1,2,...,k).
%           m =     [k x p] matrix of group means.
%

% RE Strauss, 7/13/98
%   9/7/99 -  changed sequence of output arguments.
%   1/15/00 - allow for exact reconstruction of covariances.
%   6/1/00 -  allow p=1.

function [x,g,m] = testgrps(k,n,p,exact,v,r)
  if (nargin < 1) help testgrps; return; end;

  if (nargin < 4) exact = []; end;
  if (nargin < 5) v = []; end;
  if (nargin < 6) r = []; end;

  if (isempty(v))
    v = 2;
  end;
  if (isempty(r))
    r = 0.7;
  end;
  if (isempty(exact))
    exact = 0;
  end;

  lenn = length(n);
  if (lenn~=1 & lenn~=k)
    error('  TESTGRPS: length of sample-size vector must be 1 or k');
  end;

  if (any(n<p) & exact)
    disp('  TESTGRPS warning: if number of observations is less than number of');
    disp('    variables, covariances will not be exact.');
    exact = 0;
  end;

  g = makegrps(1:k,n);              % Group identifiers

  x = [];                           % Build data matrix one group at a time
  m = zeros(k,p);

  for i = 1:k
    mu = 5 + 2*rand(1,p);
    m(i,:) = mu;

    if (lenn > 1)
      nn = n(i);
    else
      nn = n;
    end;

    if (p>1)
      x = [x; randmvn(nn,mu,v,r,exact)];
    else
      x = [x; randn(nn,1)*sqrt(v)+mu];
    end;
  end;

  return;

