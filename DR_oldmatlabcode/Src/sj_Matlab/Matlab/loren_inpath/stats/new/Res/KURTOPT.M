% KURTOPT:  Objective function for RANDKURT.  Given a vector of sample data
%           and target kurtosis value, finds the optimal value
%           of the two 'inverse' Box-Cox parameters lambda, which specify the
%           transformations (in this case) from a normal distribution to a
%           skewed and kurtotic one.  Returns the sum of squared differences
%           between the actual and desired kurtosis values.
%
%       Usage: sse = kurtopt(lambda,sampg,X)
%

% RE Strauss, 4/14/95

function sse = kurtopt(lambda,target,X)
  if (lambda < eps) | (lambda > 5)
    sse = 1e6;
    return;
  end;

  Xprime = kurttran(lambda,X);
  kurt = kurtosis(Xprime);

  if (isinf(kurt) | isnan(kurt))
    sse = 1e6;
  else
    sse = abs((target-kurt)^2);
  end;

  return;
