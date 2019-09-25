% BOXCOXF: Objective function for boxcox().  x must be strictly positive.
%         This application minimizes L rather than maximizes -L.
%

% RE Strauss, 2/10/97
%   12/8/02 - replace nu by n-1.

function L = boxcoxf(lambda,x)
  n = length(x);

  if (abs(lambda) > eps)
    xp = ((x.^lambda)-1)/lambda;
  else
    xp = log(x);
  end;

  L = -((n-1)/2)*log(var(xp)) + (lambda-1)*((n-1)/n)*sum(log(x));
  L = -L;

  return;

