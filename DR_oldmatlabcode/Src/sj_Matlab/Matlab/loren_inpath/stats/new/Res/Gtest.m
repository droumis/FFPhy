% GTEST: 2-way contingency table test using the log-likelihood G-statistic.
%
%     Usage: [G,pr,df] = gtest(obs)
%
%         obs = two-dimensional contingency table of counts.
%         --------------------------------------------------
%         G = observed value of G-statistic.
%         pr = probability of observed statistic value.
%         df = degrees of freedom for the text.
%

function [G,pr,df] = gtest(obs)
  [r,c] = size(obs);
  if (min([r,c])<2)
    error('Error: two-dimensional table needed for test');
  end;

  R = sum(obs);
  C = sum(obs')';
  N = sum(R);

  f = obs(:);
  G = 2*(sum(f.*log(f)) - sum(R.*log(R)) - sum(C.*log(C)) + N*log(N));
  df = (r-1)*(c-1);

  pr = 1-chi2cdf(G,df);

  return;
