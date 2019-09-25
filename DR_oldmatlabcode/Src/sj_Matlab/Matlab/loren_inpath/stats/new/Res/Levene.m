% LEVENE: Levene's test for unequal variability, based on anova of absolute 
%         deviations from the mean.
%
%   Syntax: [F,pr,df,ss,ms] = levene(x,grps,{iter},{CI_level})
%
%         x =         [n x 1] observations for a single variable.
%         grps =      [n x 1] classification variable.
%         iter =      optional number of randomization iterations [default = 0].
%         CI_level =  percentage confidence level for bootstrapped variance
%                       components [default=95].
%         ------------------------------------------------------------------------
%         F =       observed F-statistic value.
%         pr =      significance level of the test, either asymptotic (if iter=0)
%                     or randomized (if iter>0).
%         df =      [3 x 1] vector of degrees of freedom (among-group,
%                     within-group, total).
%         ss =      [3 x 1] vector of sums-of-squares (among-group, within-group,
%                     total).
%         ms =      [2 x 1] vector of mean-squares (among-group, within-group).
%

% RE Strauss, 5/12/99

function [F,pr,df,ss,ms] = levene(x,grps,iter,CI_level)
  if (nargin < 3) iter = []; end;
  if (nargin < 4) CI_level = []; end;

  absdev = abs(grpcentr(x,grps));     % Absolute deviations, within groups

  [F,pr,df,ss,ms] = anova(absdev,grps,iter,CI_level);

  return;

