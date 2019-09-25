% HomoMeanAIC: Estimtes homogeneous subsets of means based on Akaike 
%              information criteria, under the assumption of normality.
%              Deletes missing data.
%
%     Usage: [hs,crit,M,ordgrps] = homomeanaic(x,grps,{hetvar},{kind},{nbest})
%
%         x =       data vector for a single variable.
%         grps =    corresponding vector of group labels.
%         hetvar =  optional flag indicating whether groups variances are assumed
%                     to be homogeneous:
%                     0 = assume homogeneity of variances (single population
%                           variance) [default];
%                     1 = assume restricted heterogeneity of variances (unique
%                           population variance for each homogeneous subset).
%         kind =    optional flag indicating the information statistic to be used:
%                     0 = AIC (Akaike information criterion) [default];
%                     1 = BIC (Schwarz criterion, with penalty term of log(N));
%                     2 = CAIC (Bozdogan criterion, with penalty term of log(N)+1).
%         nbest =   optional number of best subsets to return [default = 1].
%         -------------------------------------------------------------------------
%         hs =      [nbest x k] matrix of best homogeneous subsets.
%         crit =    [nbest x 1] vector of criterion values.
%         m =       [k x 1] vector of group means.
%         ordgrps = [k x 1] corresponding vector of group identifiers.
%         

% Dayton CM. 1998. Information criteria for the paired-comparisons problem.
%   American Statistician 52:144-151.
% Dayton CM. 2001. SUBSET: best subsets using information criteria.
%   Journal of Statistical Software 6(2):1-10.
% Schwarz, G. 1978. Estimating the dimension of a model.
%   Annals of Statistics 6:461-464.
% Bozdogan, H. 1987. Model-selection and Akaike's information criterion (AIC):
%   the general theory and its analytical extensions.  Psychometrika 51:345-370.

% RE Strauss, 2/16/02

function [hs,crit,m,ordgrps] = homomeanaic(x,grps,hetvar,kind,nbest)
  if (nargin < 3) hetvar = []; end;
  if (nargin < 4) kind = []; end;
  if (nargin < 5) nbest = []; end;

  if (isempty(hetvar))
    hetvar = 0;
  end;
  if (isempty(kind))
    kind = 0;
  end;
  if (isempty(nbest))
    nbest = 1;
  end;

  if (~isvector(x))
    error('  HomoMeanAIC: data matrix must be a vector.');
  end;

  x = x(:);
  grps = grps(:);

  if (length(grps)~=length(x))
    error('  HomoMeanAIC: data and group-identification vectors not compatible.');
  end;

  i = isfinite(x);                      % Omit missing data
  x = x(i);
  grps = grps(i);

  N = length(x);                        % Total sample size
  ugrps = uniquef(grps);                % Group identifiers
  k = length(ugrps);                    % Number of groups

  bsets = dec2bin(0:((2.^k)-2));        % Binary representations of possible subsets
  bsets = bsets(1:2.^(k-1),:);
  nsets = size(bsets,1);
  gsets = ones(nsets,k);

  for i = 1:nsets                       % Convert to ordered list of partitions
    for j = 2:k
      if (bsets(i,j)~=bsets(i,j-1))
        gsets(i,j:end) = gsets(i,j:end)+1;
      end;
    end;
  end;
gsets

  nbest = min([nbest,nsets]);           % Number of best-subsets to be reported

  hs = inf*ones(nbest,1);               % Initialize output matrices
  crit = zeros(nbest,k);

  [morig,stderr,ordgrps] = means(x,grps);   % Sort means into increasing sequence
  [m,ig] = sort(morig);                 % Means
  ordgrps = ugrps(ig);                  % Corresponding group identifiers

  s2g = zeros(k,1);                     % Within-group variances (biased)
  n = zeros(k,1);                       % Sample sizes
  xx = [];                              % Sort data into corresponding sequence
  g = [];
  for i = 1:k
    j = find(grps == ordgrps(i));
    xx = [xx; x(j)];
    g = [g; grps(j)];
    n(i) = length(x(j));
    s2g(i) = var(x(j))*(n(i)-1)/n;
  end;
  x = xx;
  grps = g;

  if (hetvar==0)                        % Pooled within-group variance (biased)
    s2e = meanwt(s2g,n);
  end;

[M,stderr,V,grpids] = means(X,grps)  <-------

  for i = 1:nsets                       % Cycle thru all possible partitions
    switch(hetvar)
      case 0;                             % Homogeneous variances
        
      case 1;                             % Restricted heterogeneous variances

    end;

    ssm = 
    loglik = -(N/2)*(log(2*pi)+log(S2e)) - ssm;
  end;

  return;


