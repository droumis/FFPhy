% REGRNULL: Tests a fitted bivariate regression line against a hypothethical 
%           (null) line given by expected values of slope and intercept, taking 
%           into consideration the joint sampling distribution of b0 and b1.
%
%     Syntax: [pr,F,df,b] = regrnull(x,y,beta,{noplot})
%
%           x,y =     matching vectors (length n) for the independent and 
%                       dependent variables.
%           beta =    2-element vector of null values: [beta0, beta1].
%           noplot =  optional boolean vector indicating, if true, that a plot 
%                       of data, regression line and null line is not to be 
%                       produced [default = 0].
%           ------------------------------------------------------------------
%           pr =      significance level of test.
%           Fhat =    resulting F-statistic value.
%           df =      corresponding degrees of freedom.
%           b =       column vector of parameter estimates: [b0; b1].
%

% RE Strauss, 6/15/01
%   6/25/01 - added plot.

function [pr,Fhat,df,b] = regrnull(x,y,beta,noplot)
  if (nargin < 4) noplot = []; end;

  if (isempty(noplot))
    noplot = 0;
  end;

  x = x(:);
  y = y(:);
  beta = beta(:);

  n = length(x);
  if (length(y) ~= n)
    error('  REGRNULL: data vectors are not compatible.');
  end;
  if (length(beta) ~= 2)
    error('  REGRNULL: null-parameter vector must have two elements.');
  end;


  [b,stats] = linregr(x,y);
  mse = stats(2);
  b0 = b(1);
  b1 = b(2);
  beta0 = beta(1);
  beta1 = beta(2);

  d0 = b0 - beta0;
  d1 = b1 - beta1;
  sumx = sum(x);
  sumx2 = x'*x;

  Fhat = (n*d0*d0 + 2*sumx*d0*d1 + sumx2*d1*d1)/(2*mse);
  df = [2; n-2];
  pr = 1-fcdf(Fhat,2,n-2);

  if (~noplot)
    figure;
    plot(x,y,'ko');
    putregrline(x,y);
    hold on;
    xmin = min(x);
    xmax = max(x);
    ymin = beta0+beta1*xmin;
    ymax = beta0+beta1*xmax;
    plot([xmin,xmax],[ymin,ymax],'k--');
    hold off;
    putbnds([x;xmin;xmax],[y;ymin;ymax]);
  end;

  return;
