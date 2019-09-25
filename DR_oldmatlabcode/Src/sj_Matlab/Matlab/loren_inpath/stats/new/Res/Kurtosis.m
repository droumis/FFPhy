% KURTOSIS: Returns the unbiased estimate of kurtosis, and its standard
%           error, by column.  This measure is _not_ adjusted by subtracting 3;
%           thus the expected value for a normal distribution is 3, not zero.
%           Does not allow for missing data.
%
%     Usage: [kurt,se] = kurtosis(X)
%
%           X =    [n x p] data matrix.
%           -----------------------------------------
%           kurt = [1 x p] vector of kurtosis values.
%           se =   [1 x p] vector of standard errors.
%

% Zar 1984, pp. 81 (Eqn 7.13),119.

% RE Strauss, 4/14/95
%   6/24/99 - check for minimum sample size.

function [kurt,se] = kurtosis(X)
  [n,p] = size(X);
  if (n==1)
    X = X';
    [n,p] = size(X);
  end;

  if (n >= 4)
    sumX4 = sum(X.^4);
    sumX3 = sum(X.^3);
    sumX2 = sum(X.^2);
    sumX  = sum(X);

    t1 = (n*n*n + n*n).*sumX4;
    t2 = -4*(n*n + n).*sumX3.*sumX;
    t3 = -3*(n*n - n).*sumX2.*sumX2;
    t4 = 12*n.*sumX2.*sumX.*sumX;
    t5 = -6*(sumX.^4);
    denom = n*(n-1)*(n-2)*(n-3)*(eps+(std(X).^4));

    kurt = abs(3+((t1+t2+t3+t4+t5)./denom));

    if (nargout == 2)
      se = ones(1,p)*sqrt((24*n*(n-1)*(n-1))./((n-3)*(n-2)*(n+3)*(n+5)));
    end;
  else
    kurt = NaN*ones(1,p);
    se =   NaN*ones(1,p);
  end;

  return;
