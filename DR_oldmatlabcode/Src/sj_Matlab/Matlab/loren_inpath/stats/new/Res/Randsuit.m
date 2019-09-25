% RANDSUIT: Generate data matrix for one or more groups, allowing for character 
%           suites.
%
%     Usage: [X,suites] = randsuit(N,nsuite,{mu},{V},{wcorr},{bcorr},{exact})
%
%         N =       vector of numbers of observations per group; 
%                     sum = total sample size.
%         nsuite =  vector of numbers of characters per suite; length is the 
%                     number of suites (S), sum = number of characters (P).
%         mu =      [length N x length P] matrix of character means for each 
%                     group [default = 10 for all chars of first group, 12 for 
%                     second group, etc.].
%         V =       vector (length P) of character variances; if V is a scalar, 
%                     it is used for all characters and all groups [default = 1].
%         wcorr =   within-suite correlation (constant for all suites) 
%                     [default = 0.85].
%         bcorr =   among-suite correlation (constant for all pairs of suites) 
%                     [default = 0.65].
%         exact =   optional boolean flag indicating, if true, that covariances 
%                     among random observations are to be exactly those of the 
%                     target matrix [default = 1: covariances are exact].
%         -----------------------------------------------------------------------
%         X =       [N x P] data matrix.
%         suites =  grouping vector (length P) indicating S hypothetical 
%                     character suites.
%

% RE Strauss, 10/5/99
%   12/12/99 - added capability for multiple groups.
%   10/5/00 -  fixed problem with propagating scalar V.

function [X,suites] = randsuit(N,nsuite,mu,V,wcorr,bcorr,exact)
  if (nargin < 3) mu = []; end;
  if (nargin < 4) V = []; end;
  if (nargin < 5) wcorr = []; end;
  if (nargin < 6) bcorr = []; end;
  if (nargin < 7) exact = []; end;

  default_wcorr = 0.85;
  default_bcorr = 0.65;

  if (min(size(N))>1 | min(size(nsuite))>1)
    error('  RANDSUIT: N and nsuite must be vectors');
  end;

  nvars = sum(nsuite);
  vars = nsuite;
  nsuite = length(nsuite);
  nobs = sum(N);
  ngrps = length(N);

  if (isempty(mu))
    mu = zeros(ngrps,nvars);
    init_mu = 10;
    delta_mu = 2;
    for i = 1:ngrps
      mu(i,:) = init_mu*ones(1,nvars);
      init_mu = init_mu + delta_mu;
    end;
  end;
  if (size(mu) ~= [ngrps,nvars])
    error('RANDSUIT: mu and nsusite vectors incompatible');
  end;

  if (isscalar(V))
    V = V * ones(ngrps,nvars);
  end;
  if (isempty(V))
    V = ones(ngrps,nvars);
  end;

  if (isempty(wcorr))
    wcorr = default_wcorr;
  end;
  if (isempty(bcorr))
    bcorr = default_bcorr;
  end;
  if (isempty(exact))
    exact = 1;
  end;

  if (~exact & nobs<=nvars)
    error('  RANDSUIT: total N must be greater than number of variables');
  end;
  if (exact & min(N)<=nvars)
    disp('  RANDSUIT: for exact covariances, min N must be greater than');
    disp('              number of variables');
    disp('  Warning:  shifting to exact covariances across groups,');
    disp('              approximate covariances within groups');
    exact = 0;
  end;

  suites = makegrps([1:nsuite],[vars])';

  R = ones(nvars,nvars);
  X = zeros(sum(N),nvars);

  for i = 1:(nvars-1)
    for j = (i+1):nvars
      if (suites(i)==suites(j))
        R(i,j) = wcorr;
        R(j,i) = wcorr;
      else
        R(i,j) = bcorr;
        R(j,i) = bcorr;
      end;
    end;
  end;

  if (~exact)                             % If covars can be approximate,
    X = randmvn(nobs,zeros(1,nvars),ones(1,nvars),R,exact);   % Get random values
  else                                      % else allocate data matrix.
    X = zeros(nobs,nvars);
  end;

  lobs = 0;                               % Adjust means and variances by group
  for i = 1:ngrps
    fobs = lobs + 1;                        % First and last obs in group
    lobs = lobs + N(i);
    if (exact)
      X(fobs:lobs,:) = randmvn(N(i),mu(i,:),V(i,:),R,exact);
    else
      s = ones(N(i),1)*sqrt(V(i,:));
      m = ones(N(i),1)*mu(i,:);
      X(fobs:lobs,:) = X(fobs:lobs,:).*s + m;
    end;
  end;

  return;
