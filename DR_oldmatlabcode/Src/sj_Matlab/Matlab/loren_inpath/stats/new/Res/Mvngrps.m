% MVNGRPS:  Constructs one or more multivariate-normal samples of 
%           equal variance and covariance (within sampling variation).
%           Centroids are selected randomly from uniform distributions.
%
%           See 'randmvn' for completely specified MVN samples.
%
%     Usage: [x,g,m] = mvngrps(k,n,p,{v},{r})
%
%           k = number of groups to be created.
%           n = sample size per group (if identical for all groups), 
%                 or vector (length k) of sample sizes.
%           p = number of variables.
%           v = optional within-group variance (identical for all 
%                 variables) [default = 2].
%           r = optional within-group correlation (identical for all 
%                 pairs of variables) [default = 0.7].
%           ---------------------------------------------------------
%           g = column vector (length k*n) of group identifiers 
%                 (1,2,...,k).
%           x = [k*n x p] data matrix.
%           m = [k x p] matrix of group means.
%

% RE Strauss, 7/13/98
%   6/5/99 - rename from 'testgrps', reorder output arguments.

function [x,g,m] = mvngrps(k,n,p,v,r)
  if (nargin < 4) v = []; end;
  if (nargin < 5) r = []; end;

  if (isempty(v))
    v = 2;
  end;
  if (isempty(r))
    r = 0.7;
  end;

  lenn = length(n);
  if (lenn~=1 & lenn~=k)
    error('  MVNGRPS: length of sample-size vector must be 1 or k');
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

    x = [x; randmvn(nn,mu,v,r)];
  end;

  return;

