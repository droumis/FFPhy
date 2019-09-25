% FTEST: Given an estimated F-statistic and associated degrees of freedom,
%        returns its significance level.
%
%     Usage: prob = ftest(F,df1,df2)
%
%        F = estimated F-statistic
%        df1 = numerator degrees of freedom
%        df2 = denominator degrees of freedom

function prob = ftest(F,df1,df2)
  prob = 1-fcdf(F,df1,df2);

%  if (F<0)                         % Check for F out of range
%     error('   F-statistic must be non-negative in ftest')
%  end;

%  prob = betai(0.5*df2, 0.5*df1, df2/(df2+(df1*F)));
%  if (prob > 1)
%     prob = 1;
%  end;

  return;

