% BoxCoxInv:  Given lambda and a vector of transformed data, restores the data to its 
%             original untransformed values.  Allows for missing data.
%
%     Usage: x = boxcoxinv(xp,lambda,c)
%
%         xp =      vector of transformed values.
%         lambda =  Box-Cox parameter.
%         c =       value added to data before transforming to ensure all positive 
%                     values.
%         ----------------------------------------------------
%         x =       corresponding vector of back-transformed values.
%

% RE Strauss, 1/3/03

function x = boxcoxinv(xp,lambda,c)
  if (nargin < 2) help boxcoxinv; return; end;
  
  if (~isvector(xp))
    error('  BoxCoxInv: input data must be in vector form.');
  end;
  if (~isscalar(lambda))
    error('  BoxcoxInv: lambda must be a scalar.');
  end;
  
  i = find(isfinite(xp));
  xmean = mean(xp(i));
  xstd = std(xp(i));

  if (abs(lambda) > eps)                
    x = ((xp*lambda)+1).^(1/lambda);
  else
    x = exp(xp);
  end;
  
  x = x+c;
  
%   x = x(:);

  return;
  