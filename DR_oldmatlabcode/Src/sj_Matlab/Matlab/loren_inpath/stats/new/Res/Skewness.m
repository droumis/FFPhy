% SKEWNESS: Returns the unbiased estimate of skewness, and its standard
%           error, by column.  Doesn't allow for missing data.
%
%     Usage: [skew,se] = skewness(X)
%
%           X =    [n x p] data matrix.
%           -----------------------------------------
%           skew = [1 x p] vector of skewness values.
%           se =   [1 x p] vector of standard errors.
%

% Zar 1984, pp. 81,118.

% RE Strauss, 4/20/95
%   6/24/99 - check for minimum sample size

function [skew,se] = skewness(X)
  [n,p] = size(X);
  if (n==1)
    X = X';
    [n,p] = size(X);
  end;

  if (n >= 3)
    sum_d3 = sum((X-ones(n,1)*mean(X)).^3);
    k3 = (n*sum_d3)/((n-1)*(n-2));
    skew = k3./(eps + std(X).^3);
    skew = sign(real(skew)).*abs(skew);

    if (nargout == 2)
      se = ones(1,p)*sqrt((6*n*(n-1))/((n-2)*(n+1)*(n+3)));
    end;
  else
    skew = NaN*ones(1,p);
    se =   NaN*ones(1,p);
  end;

  return;
