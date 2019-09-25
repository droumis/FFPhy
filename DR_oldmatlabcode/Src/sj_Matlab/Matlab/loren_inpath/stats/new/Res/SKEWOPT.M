% SKEWOPT:  Objective function for RANDSTRB.  Given a vector of sample data
%           and target skewness and kurtosis values, finds the optimal value
%           of the two 'inverse' Box-Cox parameters lambda, which specify the
%           transformations (in this case) from a normal distribution to a
%           skewed and kurtotic one.  Returns the sum of squared differences
%           between the actual and target skewness and kurtosis values.
%
%       Usage: sse = skewopt(lambda,target,X)
%

% RE Strauss, 4/14/95

function sse = skewopt(lambda,target,X)
  Xprime = skewtran(lambda,X);
  skew = skewness(Xprime);

  if (isinf(skew) | isnan(skew))
    sse = 1e6;
  else
    sse = abs((target-skew)^2);
  end;
  return;
