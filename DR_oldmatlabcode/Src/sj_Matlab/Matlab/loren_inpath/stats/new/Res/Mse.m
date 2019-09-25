% MSE:  Calculate the mean-square-error (mse) and adjusted coefficient of multiple
%       determination (R^2) from observed and predicted values.  Ignores missing 
%       values.
%
%     Usage: [ms,r2] = mse(obs,pred,{p})
%
%       obs =   [n x 1] vector of observed values.
%       pred =  [n x 1] vector of corresponding predicted values.
%       p =     optional number of independent variables in the regression
%                 [default = 1].
%       ------------------------------------------------------------------
%       ms =    mean squared error.
%       r2 =    adjusted R^2.
%

% RE Strauss, 1/21/99
%   1/1/01 -  corrected degrees of freedom;
%             make p optional;
%             call r2adj() to calculate coefficient of determination.

function [ms,r2] = mse(obs,pred,p)
  if (nargin < 3) p = []; end;

  get_r2 = 0;
  if (nargout > 1)
    get_r2 = 1;
  end;

  if (isempty(p))
    p = 1;
  end;

  if (min(size(obs))~=1 | min(size(pred))~=1)
    error('  MSE: observed and predicted values must be vectors');
  end;

  n = length(obs);
  if (length(pred)~=n)
    error('  MSE: observed- and predicted-value vectors must be same length');
  end;

  i = find(~isfinite(obs) | ~isfinite(pred));   % Remove missing values
  if (~isempty(i))
    obs(i) = [];
    pred(i) = [];
    n = length(obs);
  end;

  ymean = mean(obs);
  ssto = (obs-ymean)'*(obs-ymean);    % Sums of squares
  e = obs - pred;
  sse =  e'*e;

  dfe =  n-p-1;                       % Degrees of freedom
  dfto = n-1;

  ms =  sse / dfe;
  msto = ssto / dfto;

  if (get_r2)
    ssr = ssto - sse;
    r2 = r2adj(ssr,ssto,n,p);         % Adjusted coeff of multiple determination
  end;

  return;

