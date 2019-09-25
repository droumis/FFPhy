% ANOVAF: Objective function for ANOVA, called by BOOTSTRP.
%
%     Syntax: solution = anovaf(x,grps,nu1,nu2,getF,getcomps)
%
%         x =        [n x p] data matrix.
%         grps =     [n x 1] vector of group memberships.
%         nu1,nu2 =  unused variables provided by BOOTSTRP.
%         getF =     boolean flag indicating that F-value is to be returned.
%         getcomps = boolean flag indicating that variance components are to
%                      be returned.
%                    NOTE: Only one of the two flags can be true.
%         ------------------------------------------------------------------
%         solution = F statistic (if getF==1), or
%                    [1 x 4] vector of [varcomp varprop] (if getcomps==1).
%

% RE Strauss, 9/26/96
%   1/16/00 - minor changes to documentation and error message.

function solution = anovaf(x,grps,nu1,nu2,getF,getcomps)

  n = length(x);                      % Number of observations
  G = design(grps);                   % Design matrix
  ngrps = size(G,2);                  % Number of groups
  totmean = mean(x)*ones(n,1);        % Matching vector of grand means
  y = x - totmean;                    % Deviations from grand mean
  ssto = y'*y;                        % Total sum-of-squares

  dfto = n-1;                         % Total df
  dfa  = ngrps-1;                     % Among-group df
  dfe  = dfto - dfa;                  % Within-group df

  gpg = G'*G;
  gbar = inv(gpg)*G'*x;               % Group means
  xbar = G*gbar;                      % Matching vector of group means
  e = x - xbar;                       % Deviations from group means

  sse = e'*e;                         % Within-group sum-of-squares
  ssa = ssto - sse;                   % Among-group sum-of-squares

  msa = ssa / dfa;                    % Among-group mean-squares
  mse = sse / dfe;                    % Within-group mean-squares

  % F-statistic

  if (getF)
    solution = msa / mse;               % Observed F-statistic
  end;

  % Estimates of variance components

  if (getcomps)
    gn = diag(gpg);                     % Group sample sizes
    n0 = (1/(ngrps-1)) * (sum(gn)-(gn'*gn/sum(gn)));
    s2e = mse;
    s2a = (msa-mse)/n0;
    s2to = s2a + s2e;
    varcomp = [s2a, s2e];               % Variance components
    varprop = [s2a/s2to, s2e/s2to];     % Among & within proportions
    solution = [varcomp varprop];
  end;

  return;


