% TobitLLFn: log-likelihood function for Tobit regression.
%
%     Usage: L = tobitllfn(p,X,y)
%
%         p = vector of parameter values: sigma and beta estimates.
%         X = [n x p] matrix of independent variables.
%         y = [n x 1] vector of dependent variables, with NaNs for censored values.
%         ---------------------------------------------------------------------------
%         L = -(log-likelihood), to be minimized for the maximum-likelihood solution.
%

% RE Strauss, 10/28/02
%   8/11/03 - protect against s==0.
%   8/14/03 - protect against phi==0 or phi==1.

% Breen, R. 1996. Regression models: censored, sample selected, or truncated data.
%   Sage Publs.

function L = tobitllfn(p,X,y)
  s = max([abs(p(1)),eps]);             % Extract parameters
  s2 = s*s;
  b = p(2:end)';
  
  Xb = X*b;                             % Predicted values
  
  i = find(~isfinite(y));               % Censored observations
  phi = normcdf(Xb(i,:)/s);
  i = find(phi==0);
  if (~isempty(i))
    phi(i) = eps;
  end;
  i = find(phi==1);
  if (~isempty(i))
    phi(i) = 1-eps;
  end;
  lnL = sum(log(1-phi));                
  
  i = find(isfinite(y));                % Uncensored observations
  lnL = lnL + length(i)*log(1/sqrt(2*pi*s2));             
  lnL = lnL - sum((1/(2*s2))*(y(i)-Xb(i,:)).^2);

  L = -lnL;                             % Minimize rather than maximize
  return;
  