% Chi2Comb: Given a set of chi-squared values for a single value of degrees-of-freedom,
%           combines them into an F-statistic to provide a single p-value.
%
%     Usage:  [p,F,df1,df2] = chi2comb(chi2,df)
%
%         chi2 =    vector of chi-squared values.
%         df =      scalar degrees-of-freedom.
%         -----------------------------------------------
%         p =       overall significance level for tests.
%         F =       F-statistic value.
%         df1,df2 = degrees of freedom for F-test.
%

% Schafer, J.L. (1997) Analysis of Incomplete Multivariate Data.  London: 
%   Chapman and Hall, p. 15.

% RE Strauss

function [p,F,df1,df2] = chi2comb(chi2,df)
  chi2 = chi2(:);
  if (~isscalar(df))
    error('  CHI2COMBINE: df must be a scalar.');
  end;

  m = length(chi2);
  g = sqrt(chi2);
  mchi2 = mean(chi2);
  r = (1+1/m)*((g'*g)-(sum(g).^2)/m)/(m-1);
  F = (mchi2/df - r*(m-1)/(m+1))/(1+r);
  df1 = df;
  df2 = (m-1)*(1+1/r)^2/df^(3/m);
  p = 1 - fcdf(F,df1,df2);
  
  return;
 