% TStatistic: Calculates the value of the t-statistic for two groups.
%
%     Usage: t = tstatistic(x,g)
%
%         x = [n x 1] vector of data for a single variable.
%         g = corresponding group-membership vector.
%         -------------------------------------------------
%         t = value of the t-statistic.
%

function t = tstatistic(x,g)
  [grpid,n] = uniquef(g);
  
  x1 = x(g==grpid(1));
  x2 = x(g==grpid(2));
  
  m1 = mean(x1);
  m2 = mean(x2);
  
  v1 = var(x1);
  v2 = var(x2);
  
  t = (m1-m2)./(sqrt(v1/n(1) + v2/n(2)));

  return;
  