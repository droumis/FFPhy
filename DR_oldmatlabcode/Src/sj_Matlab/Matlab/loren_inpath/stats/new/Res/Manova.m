% MANOVA: Single-classification multivariate ANOVA for k groups and p variables.
%         Uses the F approximation to Wilks' lambda.
%         Optionally randomizes the observations among groups to determine
%         the significance level of lambda.
%
%   Syntax: [lambda,F,pr,df,Fdf] = manova(X,grps,{iter},{CI_level})
%
%         X =         [n x p] analytic variable.
%         grps =      [n x 1] classification variable.
%         iter =      optional number of randomization iterations to estimate
%                       significance level (default = 0).
%         CI_level =  percentage confidence level for bootstrapped variance
%                       components (default = 95).
%         ----------------------------------------------------------------------
%         lambda =    observed Wilks' lambda.
%         F =         estimated F-statistic value.
%         pr =        significance level of the test, either asymptotic
%                       (if iter=0) or randomized (if iter>0).
%         df =        [1 x 3] vector of degrees of freedom (among-group,
%                       within-group, total).
%         Fdf =       [1 x 2] vector of degrees of freedom for the F test (numerator, 
%                       denominator).
%

% Tabachnick, BG & LS Fidell. 1989. Using Multivariate Statistics (2nd ed.),
%   pp. 381-389.  Harper-Collins.
% Seber, GAF.  1984.  Multivariate Observations, pp. 433-437.  Wiley Series in 
%   Probability and Mathematical Statistics.

% RE Strauss
%   11/9/98 -   error message for singular E+H matrix.
%   11/29/99 -  changed calling sequence.
%   6/13/00 -   added check for missing data.

function [lambda,F,pr,df,Fdf] = manova(X,grps,iter,CI_level)
  if (nargin < 3) iter = []; end;
  if (nargin < 4) CI_level = []; end;

  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(CI_level))
    CI_level = 0.95;
  elseif (CI_level > 1)
    CI_level = CI_level/100;
  end;
  alpha = 1-CI_level;

  [N,nvars] = size(X);
  p = nvars;

  if (misscheck(X,grps))
    error('  MANOVA: data matrix or grouping vector contains missing data.');
  end;

  G = design(grps);                   % Design matrix
  ngrps = size(G,2);                  % Number of groups

  totmean = ones(N,1)*mean(X);        % Matching matrix of grand means

  mean_W = (G'*G)\G'*X;               % Within-group means
  grpmean = G*mean_W;                 % Matching matrix of group means

  e = X - grpmean;                    % Within-group deviations
  g = totmean - grpmean;              % Among-group deviations
% y = X - totmean;                    % Total deviations

  sse =  e'*e;                        % Within-group sscp
  ssa =  g'*g;                        % Among-group sscp
% ssto = y'*y;                        % Total sscp

  dfto = N-1;                         % Total df
  dfa  = ngrps-1;                     % Among-group df
  dfe  = dfto - dfa;                  % Within-group df
  df = [dfa dfe dfto];

  H = ssa;                            % Hypothesis matrix
  E = sse;                            % Error matrix

  rankEH = rank(E+H);
  sizeH = size(H,1);
  if (rankEH < size(H,1))
    disp(sprintf('  MANOVA warning: E+H matrix is singular (rank %1.0f < %1.0f)',...
                    rankEH,sizeH));
  end;

  lambda = det(E)/(det(E+H));         % Wilks' lambda

  if (p==2 & ngrps==2)                % Fix for single case in which fails
    s = 2;
  else
    denom = ((p*p)+(dfa*dfa)-5);
    if (denom>0)
      s = sqrt(((p*p*dfa*dfa)-4)/denom);
    else
      s = NaN;
    end;
  end;

  y = lambda.^(1/s);

  df1 = floor(p*dfa);
  df2 = floor(s*(dfe-((p-dfa+1)/2))-((p*dfa-2)/2));
  Fdf = [df1 df2];

  F = ((1-y)/y)*(df2/df1);

  % Significance level of observed F-statistic

  if (iter == 0)                      % Asymptotic significance level
    pr = 1 - fcdf(F,df1,df2);
  else                                % Randomized significance level
    [ci,pr] = bootstrp('manovaf',[0,1,0],iter,alpha,X,grps,0);
  end;

  return;
