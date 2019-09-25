% MAHALPR: Asymptotic probabilities for Mahalanobis distances among groups,
%          based on a central F-distribution.
%
%     Syntax: prob = mahalpr(X,grps,D2)
%
%        X =         [n x p] data matrix (obs x vars).
%        grps =      row or column vector of group identifiers.
%        D2 =        [g x g] matrix of Mahalanobis distances.
%        ---------------------------------------------------------------------
%        prob =      [g x g] matrix of probabilities of H0:D2=0.
%

% SPSS 7.5 Statistical Algorithms, p 136.

% RE Strauss, 5/6/98
%   11/29/99 - changed calling sequence.

function prob = mahalpr(X,grps,D2)
  [index,N] = uniquef(grps);
  g = length(index);            % Number of groups
  [nobs,q] = size(X);           % Numbers of observations & variables

  n = N*ones(1,g) + ones(g,1)*N';   % Sample size for each pairwise comparison
  c = (N*N')./n;
  m = n-q-g+1;
  k = q*(n-q);  

  df1 = q;
  df2 = n-q-1;

  Fhat = D2.*m.*c./k;
  prob = 1 - fcdf(Fhat,df1,df2);

  return;

