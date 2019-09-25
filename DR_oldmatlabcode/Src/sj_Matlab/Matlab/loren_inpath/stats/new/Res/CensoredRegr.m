% CensoredRegr: Fits the regression model for censored data using the Tobit model, in 
%         which all values of the predictor variable(s) are known but the response variable 
%         is truncated below a fixed censor point,.  Assumes that censored values are 
%         passed as NaN.  Calls 'tobitregr' to fit the model, but bootstraps the standard
%         errors and covariances of the regression coefficients.
%
%     Usage: [b,stats,bstderr,bcorr,ynew,ypred] = censoredregr(X,y,{limit},{iter},{doplot})
%
%         X =       [n x p] vector of independent-variable (predictor) values.
%         y =       corresponding vector of dependent-variable (response) values, with 
%                     censored observations passed as NaN.
%         limit =   optional censoring limit value for dependent variable [default = minimum
%                     non-censored value - eps].
%         iter =    optional number of bootstrap iterations [default = 0].
%         noplot =  optional boolean variable indicating, if true, that a scatterplot is
%                     to be produced of y vs a predictor variable.  If doplot==0, no plot
%                     is produced.  If doplot>0, a plot of y vs X(:,doplot) is produced;
%                     i.e., the value of doplot indicates which predictor variable is
%                     plotted.  [Default = 1].
%         ---------------------------------------------------------------------------------
%         b =       regression coefficients [b0, b1, ..., bp]'.
%         stats =   vector of regression statistics based only on non-censored data:
%                     1: adjusted coefficient of determination (R_a^2);
%                     2: MSE (estimated residual variance);
%                     3: F-statistic value;
%                     4: df1 (numerator degrees of freedom, =dfr);
%                     5: df2 (denominator degrees of freedom =dfe);
%                     6: pr >= F (probability under the null hypothesis).
%         bstderr = corresponding vector of bootstrapped standard errors of b.
%         bcorr =   bootstrapped sampling correlations among coefficients.
%         ynew =    corresponding vector of actual response values for non-censored points
%                     and predicted response values for censored points.
%         ypred =    predicted values, including predicted censored values.
%

% RE Strauss, 12/4/02
%   12/5/02 - added optional plot.
%   12/6/02 - refined the plot; added vector of regression statistics.

function [b,stats,bstderr,bcorr,ynew,ypred] = censoredregr(X,y,limit,iter,doplot)
  if (~nargin) help censoredregr; return; end;

  if (nargin < 3) limit = []; end;
  if (nargin < 4) iter = []; end;
  if (nargin < 5) doplot = []; end;
  
  if (nargout < 3) iter = 0; end;
  
  if (isempty(iter)) iter = 0; end;
  
  [n,p] = size(X);
  [m,q] = size(y);
  if (~isvector(y))
    error('  CensoredRegr: dependent variable must be a vector.');
  end;
  if (n~=m)
    error('  CensoredRegr: X and y matrices must have same sample size.');
  end;
  if (limit > max(y))
    error('  CensoredRegr: data must be censored from below.');
  end;
  if (limit > min(y))
    error('  CensoredRegr: dependent-variable values must be > censoring limit.');
  end;
  
  if (isempty(limit))                       % Upper estimate on censoring limit
    limit = min(y(isfinite(y))) - eps;
  end;
  if (isempty(doplot))
    doplot = 1;
  end;

  ic = find(~isfinite(y));                  % Indices of censored values
  inc = find(isfinite(y));                  % Indices of non-censored values
  nnc = length(inc);                        % Number of non-censored values
  
  [b,mse,ypred] = tobitregr(X,y,limit);     % Fitted regression model
  p = length(b)-1;
  
  bstderr = [];
  bcorr = [];
  
  if (iter>0)
    bresults = zeros(iter,p+1);             % Bootstrap coefficients
    for it = 1:iter
      byX = bootsamp([y X]);
      by = byX(:,1);
      bX = byX(:,2:end);
      if (any(~isfinite(by)))
        bb = tobitregr(bX,by,limit);
      else
        bb = linregr(bX,by);
      end;
      bresults(it,:) = bb';
    end;
    
    bstderr = std(bresults)';
    bcorr = corrcoef(bresults);
  end;
  
  dfr = p;                                  % Degrees of freedom
  dfe = nnc-p-1;

  % Note: the following inflation of SSTO gives the same R2a as the
  % regression of the non-censored data.
  ssto = var(y(inc))*(nnc-1);               % Total ssq for non-censored data
  varcensX = sum(diag(var(X(inc,:))));
  varallX = sum(diag(var(X)));
  ssto = ssto * (varallX / varcensX);       % Inflate to account for full range of X
  
  e = y(inc)-ypred(inc);                    % Residuals for non-censored data
  sse =  e'*e;
  ssr =  ssto - sse;                        % Regression ssq for non-censored data  

  [r2a,F,prF] = r2adj(ssr,ssto,nnc,p+1);    % Adjusted R^2 and F-test
  stats = [r2a mse F dfr dfe prF]';         % Accumulate stats
  
  ynew = y;
  ynew(ic) = ypred(ic);   
    
  if (doplot)
    x = X(:,doplot);
    [xmin,imin] = min(x);
    [xmax,imax] = max(x);
    ypredmin = ypred(imin);
    ypredmax = ypred(imax);
    [ub,ustats,upred] = linregr(x(inc),y(inc),[],[xmin,xmax]);
    
    plot(...
      [xmin,xmax],upred,'k:',...
      [xmin,xmax],[ypredmin,ypredmax],'k',...
      x,ynew,'ko',...
      x(ic),ynew(ic),'k*',...
      [xmin,xmax],[limit,limit],'k--');
    putbnd(x,ynew);
    legend('Uncensored only','Tobit','Uncensored','Predicted censored',2);
    putxlab(sprintf('Predictor %d',doplot));
    putylab('Response');
  end;
  
  return;
  