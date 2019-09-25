% MAHALCI: Asymptotic confidence interval for Mahalanobis distances among groups, 
%          based on non-central F distribution.
%
%     Syntax: CI = mahalci(X,grps,D2,{CI_level})
%
%        X =         [n x p] data matrix (obs x vars).
%        grps =      row or column vector of group identifiers.
%        D2 =        [g x g] matrix of Mahalanobis distances.
%        CI_level =  percent width of confidence intervals [default=95].
%        ---------------------------------------------------------------------
%        CI =    [g x g] matrix of low (lower triangular matrix) and
%                  high (upper triangular matrix) confidence limits.
%                  If iter=0, they are asymptotic estimates using a
%                    noncentral F-distribution.
%                  If iter>0, they are unbiased bootstrapped estimates.
%

% RE Strauss, 5/6/98
%   11/29/99 - changed calling sequence.

function CI = mahalci(X,grps,D2,CI_level)
  index = uniquef(grps);
  ngrps = length(index);        % Number of groups
  [nobs,p] = size(X);           % Numbers of observations & variables
  ndists = ngrps*(ngrps-1)/2;   % Pairwise combinations

  if (nargin < 4)
    CI_level = [];
  end;

  if (isempty(CI_level))        % If CI-level not passed, set to default
    CI_level = 0.95;
  else
    if (CI_level > 1)             % Convert to interval [0,1]
      CI_level = 0.01*CI_level;
    end;
  end;

  CI = zeros(ngrps,ngrps);      % Confidence intervals
  N =  zeros(ngrps,1);          % Within-group sample sizes

  G = design(grps);             % Design matrix
  mean_W = (G'*G)\G'*X;         % Within-group means
  for g = 1:ngrps               % Pooled within-group covariance matrix
    index = find(G(:,g));         % Indices to nonzero elements
    N(g) = length(index);         % Sample size of current group
  end;

  sumN = N*ones(1,ngrps) + ones(ngrps,1)*N';
  k = sumN-2;
  m = k-p+1;
  c = (N*N') ./ sumN;

  sumrecipN = (1./N)*ones(1,ngrps) + ones(ngrps,1)*(1./N)';

  delta2 = (m.*D2./k)-p*(sumrecipN);    % Anderson 1984:219
%    lambda = c .* delta2 ./ k;             % Krzanowski & Marriott
  lambda = c .* delta2;                  % Muirhead
%    lambda = delta2;                  % My derivation

%Fmean = m.*(p+lambda)./((m-2).*p)
%lambda_mean = Fmean.*(m-2).*p./(m.*c) - p./c

  Fhat = zeros(ngrps,ngrps);

  limlow = (1-CI_level)/2;
  limhigh = 1 - limlow;

  for g1 = 1:(ngrps-1)
    for g2 = (g1+1):ngrps
      m12 = m(g1,g2);
      if (m12 > 1)
%          k12 = k(g1,g2);
%          c12 = c(g1,g2);
        lambda12 = lambda(g1,g2);
%          d12 = D2(g1,g2);

%pr = betanc(x,p/2,m12/2,lambda12);

        Fhat(g1,g2) = finvnc(limhigh,p,m12,lambda12);
        Fhat(g2,g1) = finvnc(limlow,p,m12,lambda12);
      end;
    end;
  end;

  FtoD2 = p.*k./(m.*c);
  CI = Fhat .* FtoD2;

  return;

