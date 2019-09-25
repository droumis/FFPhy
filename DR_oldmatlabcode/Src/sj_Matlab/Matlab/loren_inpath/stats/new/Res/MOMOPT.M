% MOMOPT:  Objective function for RANDSTRB.  Given a vector of sample data
%           and target skewness and kurtosis values, and values of
%           the two 'inverse' Box-Cox parameters lambda that specify the
%           transformations (in this case) from a normal distribution to a
%           skewed and kurtotic one, returns the sum of squared differences
%           between the actual and desired skewness and kurtosis values.
%
%       Usage: sse = momopt(lambda,target,X)
%

% RE Strauss, 5/26/95

function sse = momopt(lambda,target,X)
  if (lambda(2) <= 0) | (lambda(2) > 5)
    sse = 1e6;
    return;
  end;

%  X = skewtran(lambda(1),X);
  X = kurttran(lambda(2),X);
  X = skewtran(lambda(1),X);
  s = std(X);                         % Transformed standard deviation
  X = (target(3)/s)*X;                % Final rescale

  skew = skewness(X);
  kurt = kurtosis(X);

  if (isinf(skew) | isnan(skew) | ~isreal(skew) | ...
      isinf(kurt) | isnan(kurt) | ~isreal(kurt))
    sse = 1e6;
    return;
  end;

  sse = abs((target(1)-skew)^2 + (target(2)-kurt)^2);

%disp([lambda(1) lambda(2) skew kurt sse]);
  return;
