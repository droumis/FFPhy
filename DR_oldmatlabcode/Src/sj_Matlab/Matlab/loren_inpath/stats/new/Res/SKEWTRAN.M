% SKEWTRAN:
%
%       Usage: Xprime = skewtran(lambda,X)
%

% Box, Hunter & Hunter (1978:239)

% RE Strauss, 5/26/95

function Xprime = skewtran(lambda,X)
  gm = exp(mean(log(X)));             % Geometric mean

  if (abs(lambda) < eps)
    Xprime = gm * log(X);
  else
%    Xprime = sign(lambda)*((X.^lambda)/(abs(mean(X)).^lambda));
    Xprime = ((X.^lambda)-1)./(lambda * gm^(lambda-1));
  end;
  return;
