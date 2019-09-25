% CanCorrRegr: Canonical correlation regression of a suite of correlated independent variables
%              on a single dependent variable.  
%
%     Usage: [b,scoef,ccoef,R2a,mse,pr,df] = cancorrregr(X,y)
%
%         X =     [n x p] matrix of independent variables.
%         y =     [n x 1] vector of independent variable.
%         -------------------------------------------------------------------------------
%         b =     vector of coefficients ([b0 b1]') from regression of y on cancorr scores.
%         scoef = scoring coefficients of independent variables on canonical variate.
%         ccoef = vector correlations of independent variables with canonical variate.
%         R2a =   [1 x 2] vector of adjusted R^2 from regression of y on cancorr scores:
%                   1: original adjusted R^2 from regression of y on cancorr scores;
%                   2: adjusted R^2 from jackknifed residuals.
%         mse =   [1 x 2] vector of mean-square-errors from regression of y on cancorr scores
%                 ( original and jackknifed).
%         pr =    [1 x 2] vector of corresponding significance levels (original and jackknifed).
%         F =     [1 x 2] vector of corresponding F statistics (original and jackknifed).
%         df =    [1 x 2] vector of numerator and denominator degrees of freedom.
%

% RE Strauss, 11/8/02
%   11/14/02 - added pr,F,df as output arguments.

function [b,scoef,ccoef,R2a,mse,pr,F,df] = cancorrregr(X,y)
  y = y(:);
  n = length(y);
  if (size(X,1)~=n)
    error('  CanCorrRegr: input matrices X & y not compatible.');
  end;
  
  loadtype = 2;                       % Scoring coefficients
  [r,scoef,y_load,x_scores,y_scores] = cancorr(X,y,[],loadtype);
  [r,ccoef] = cancorr(X,y);           % Correlations of variables with canonical variate
  scoef = scoef(:,1);
  ccoef = ccoef(:,1);
  [b,stats] = linregr(X,y_scores);
  
  resid = zeros(n,1);
  for ij = 1:n                        % Jackknife the R2a and mse statistics
    i = 1:n;
    i(ij) = [];
    Xj = X(i,:);
    yj = y(i);
    [r,x_load,y_load,x_scores,y_scores] = cancorr(Xj,yj,[],loadtype);
    [bj,statsj,predj,resid(ij)] = linregr(Xj,y_scores,[],X(ij,:),y(ij));
  end;
  
  ssto = var(y)*(n-1);
  sse = var(resid)*(n-1);
  ssr = ssto-sse;
  p = length(b);
  [jr2adj,jF,jpr,df] = r2adj(ssr,ssto,n,p);
  
  R2a = [stats(1), jr2adj];
  mse = [stats(2), var(resid)];
  pr =  [stats(6), jpr];
  F =   [stats(3), jF];
  
  return;
  