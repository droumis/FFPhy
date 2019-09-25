% TobitRegr: Fits the regression model for censored data, in which all values of the
%         independent variable are known but the dependent variable is truncated below a fixed
%         censor point.  Assumes that censored values are passed as NaN.
%
%     Usage: [b,mse,pred] = tobitregr(X,y,{limit})
%
%         X =     [n x p] vector of independent-variable (predictor) values.
%         y =     corresponding vector of dependent-variable (response) values, with 
%                   censored observations passed as NaN.
%         limit = optional censoring limit value for dependent variable [default = minimum
%                   non-censored value - eps].
%         ---------------------------------------------------------------------------------
%         b =     regression coefficients [b0, b1, ..., bp].
%         mse =   mean square error (estimated residual variance).
%         pred =  predicted values, including predicted censored values.
%

% Note: add adjusted R-squred, p-value.

% RE Strauss, 10/28/02
%   12/4/02 - changed default value of 'limit' from zero to min non-censored value of y
%               minus a small value.
%   12/5/02 - added 'limit' to predicted values.

function [b,mse,pred] = tobitregr(X,y,limit)
  if (nargin < 3) limit = []; end;
  
  [n,p] = size(X);
  [m,q] = size(y);
  if (~isvector(y))
    error('  TobitRegr: dependent variable must be a vector.');
  end;
  if (n~=m)
    error('  TobitRegr: X and y matrices must have same sample size.');
  end;
  if (limit > max(y))
    error('  TobitRegr: data must be censored from below.');
  end;
  if (limit > min(y))
    error('  TobitRegr: dependent-variable values must be > censoring limit.');
  end;
  
  if (isempty(limit))
    limit = min(y(isfinite(y))) - eps;
  end;

  b = [];
  
  X = [ones(n,1),X];                % Pre-concatenate a vector of one's for intercept estimate
  y = y-limit;
  
  i = find(isfinite(y));            % Uncensored observations
  Xu = X(i,:);
  invX = inv(Xu'*Xu);
  b = invX*Xu'*y(i);                % Initial regression estimates

  pred = X*b;
  e = y(i) - pred(i);                         
  dfe = n-p-1;
  if (dfe<1)
    error('  TobitRegr: error degress of freedom are <1');
    return;
  end;
  mse =  (e'*e)/dfe;
  sigma = sqrt(mse);                % Initial estimate of sigma

  p = fminsearch('tobitllfn',[sigma,b'],[],X,y);

  sigma = p(1);
  b = p(2:end)';
  
  pred = X*b+limit;                 % Final estimates
  b(1) = b(1)+limit;                % Add censor limit to intercept
  mse = sigma*sigma;
  
  return;
  