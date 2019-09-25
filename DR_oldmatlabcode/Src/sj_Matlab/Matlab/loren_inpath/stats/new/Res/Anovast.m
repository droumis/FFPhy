% ANOVAST: Performs a 1-way anova on estimates of a statistic S, given 
%          the standard errors of the estimates and the samples sizes 
%          that were used in their estimation.
%
%     Syntax: [F,pr,df,ss,ms] = anovast(S,stderr,n)
%
%          S =       vector (length k) of statistic estimates for k groups.
%          stderr =  corresponding vector of standard errors.
%          n =       corresponding vector of sample sizes.
%          ----------------------------------------------------------
%          F =       observed F-statistic value.
%          pr =      significance level of the test.
%          df =      [1 x 3] vector of degrees of freedom (among-group,
%                      within-group, total).
%          ss =      [1 x 3] vector of sums-of-squares (among-group, within-group,
%                      total).
%          ms =      [1 x 2] vector of mean-squares (among-group, within-group).
%

% RE Strauss, 2/18/97

function [F,pr,df,ss,ms] = anovast(S,stderr,n)

  if (min(size(S))>1 | min(size(stderr))>1 | min(size(n))>1)
    error('  ANOVAST: input arguments must be vectors.');
  end;

  k = length(S);
  if (length(stderr)~=k | length(n)~=k)
    error('  ANOVAST: input vectors must be of same length.');
  end;

  x = zeros(sum(n),1);
  grps = x;

  first = 1;
  for g = 1:k                   % Create samples for groups
    ng = n(g);
    last = first + ng-1;
    n1 = ng+1;
    perc = [(1/n1):(1/n1):1-(1/n1)]';
    x(first:last) = S(g) + tinv(perc,ng-1)*stderr(g)*sqrt(ng);
    grps(first:last) = g*ones(ng,1);
    first = last+1;
  end;

%plot(grps,x,'o');
%putbnd(grps,x);

  [F,pr,df,ss,ms] = anova(x,grps);

  return;
