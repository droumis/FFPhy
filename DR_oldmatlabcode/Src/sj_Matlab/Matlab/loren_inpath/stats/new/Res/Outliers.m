% Outlier detection by measuring the difference between observed and expected

% Problem: the jpred are systematically too large after the inverse Box-Cox transform.  Why?


function outliers(x)
  x = x(:);
  x = sort(x);
  i = find(~isfinite(x));
  x(i) = [];

  [n,p] = size(x);
  
  [xp,lambda,c] = boxcoxnorm(x,0);
  z = rankits(n);
  
  [b,stats,pred,resid,stderrs,B,jpred] = linregr(z,xp);
  jpred = boxcoxinv(jpred,lambda,c);
  
  mn = min([max(x),max(jpred)]);
  mx = max([min(x),min(jpred)]);
 
  scatter(x,jpred);
  putxlab('Observed data');
  putylab('Jackknifed predicted values');
%   hold on;
%   plot([mn mx],[mn mx],'k');
%   hold off;
putregrline(x,jpred);

  return;
  