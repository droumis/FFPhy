% LINREGRF: Linear regression calculations for linregr().
%
%     Usage: [b,sb,B,stats,e] = linregrf(X,y,W,weighted)
%

% RE Strauss, 12/18/01 - extracted from linregr().

% 9/29/02 - added standard errors of estimates and beta (standardized) coefficients.
% 10/1/02 - remove conditional computation of regression statistics;
%           remove 'calc_stats' input argument.
% 1/1/03 -  guard against division by zero in calculation of beta coefficients.

function [b,sb,B,stats,e] = linregrf(X,y,W,weighted,calc_stats)
  [n,p] = size(X);
  p = p-1;                            % Omit constant from number of parameters

  if (weighted)                       % Multiple regression
    invX = inv(X'*W*X);
    b = invX*X'*W*y;      
  else
    invX = inv(X'*X);
    b = invX*X'*y;
  end;
  
  pr = X*b;                           % Predicted values
  e = y - pr;                         % Residuals
  B = NaN*ones(size(b));              % Prep for standardized coefficients

  dfr = p;                            % Degrees of freedom
  dfe = n-p-1;
  dfto = n-1;

  if (dfe<1)
    dfe = NaN;
  end;

  ssto = var(y)*dfto;                 % Sums of squares
  sse =  e'*e;
  ssr =  ssto - sse;

  msr =  ssr / dfr;                   % Mean squares
  mse =  sse / dfe;
  msto = ssto / dfto;
    
  sb = sqrt(diag(invX)*mse);          % Stderrs of estimates
    
  varX = var(X)';                     % Standardized coefficients
  vary = var(y)*ones(p+1,1);
  if (all(vary>eps))
    B = b.*sqrt(varX./vary);
  end;

  if (any(p+1>=n) | abs(ssto-ssr)<eps)
    r2a = NaN;
    F = NaN;
    prF = NaN;
  else
    [r2a,F,prF] = r2adj(ssr,ssto,n,p+1);   % Adjusted R^2 and F-test
  end;

  stats = [r2a mse F dfr dfe prF]'; % Accumulate stats

  return;
  