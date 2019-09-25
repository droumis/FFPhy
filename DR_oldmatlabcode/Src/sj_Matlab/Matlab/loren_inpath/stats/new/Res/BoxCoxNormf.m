% BoxCoxNormf: Objective function for boxcoxnorm().  x must be strictly positive.
%         This application minimizes -W rather than maximizes W.

% RE Strauss, 1/1/03, modified from boxcoxf().

function W = boxcoxf(lambda,x,censdir)
  if (abs(lambda) > eps)
    xp = ((x.^lambda)-1)/lambda;
  else
    xp = log(x);
  end;

  W = normaltest(real(xp),0,censdir);
  W = -W;
  
  return;

