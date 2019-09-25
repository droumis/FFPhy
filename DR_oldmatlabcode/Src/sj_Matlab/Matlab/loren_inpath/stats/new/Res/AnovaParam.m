%  AnovaParam: ANOVA for k groups, given only the parameters (means, standard 
%       deviations, and sample sizes) for the k groups.
%
%     Usage: [F,df,pr] = anovaparam(m,s,n)
%
%         m =  [k x 1] vector of group means.
%         s =  corresponding vector of standard deviations.
%         n =  corresponding vector of sample sizes.
%         --------------------------------------------------------------------------
%         F =  F-statistic value.
%         df = [3 x 1] vector of degrees of freedom (numerator, denominator, total).
%         pr = significance level.
%

% RE Strauss, 10/15/02

function [F,df,pr] = anovaparam(m,s,n)
  m = m(:);
  s = s(:);
  n = n(:);

  k = length(m);
  N = sum(n);
  
  if (length(s)~=k | length(n)~=k)
    error('  AnovaParam: input vectors do not match.');
  end;

  grandmean = meanwt(m,n);
  sm2 = sum(n.*(m-grandmean).^2)/(k-1);
  ssd = (s.*s).*(n-1);
  sp2 = sum(ssd)/(N-k);
  
  F = sm2/sp2;
  df = [k-1; N-k; N-1];
  pr = 1 - fcdf(F,df(1),df(2));

  return;
  