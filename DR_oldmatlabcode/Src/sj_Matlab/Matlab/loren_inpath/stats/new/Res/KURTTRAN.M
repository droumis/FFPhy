% KURTRAN:  Given a distribution vector X and a 'modified' Box-Cox parameter
%           lambda, alters the kurtosis of the distribution by transforming the
%           deviations from the mean.
%
%       Usage: Xprime = kurttran(lambda,X)
%

% RE Strauss, 5/26/95

function Xprime = kurttran(lambda,X)
  m = mean(X);
  d = X-m;
  ad = abs(d);

  if (abs(lambda) < eps)
%    Xprime = m + sign(d).*log(ad);
%    lambda = eps;

    Xprime = m + sign(d).*log(ad+1);
  else
%    Xprime = exp(lambda*ad/mean(ad));
%    Xprime = m + sign(d).*(Xprime-min(Xprime));
%    Xprime = m + sign(d).*((ad.^lambda)./lambda);

    Xprime = m + sign(d).*((((ad+1).^lambda)-1)./lambda); % John & Draper (1980)
  end;

  return;
