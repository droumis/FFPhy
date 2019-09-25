% CHI2NC: Evaluates the cumulative probability of a noncentral chi-squared
%         distribution with parameter df and noncentrality lambda.
%
%     Syntax: prob = chi2nc(chi2,df,lambda)
%

function prob = chi2nc(chi2,df,lambda)
  dl = df + lambda;
  d2l = df + 2.*lambda;
  dl2 = 9.*dl.*dl;

  z = (chi2./dl).^(1/3) - 1;
  z = z + 2.*d2l./dl2;
  z = z ./ sqrt(2.*d2l./dl2);

  prob = normcdf(z);

  return;

