% MEANWT: Calculate a weighted mean and (population) variance, skewness, and 
%         kurtosis values.
%
%     Usage: [wt_mean,wt_var,wt_skew,wt_kurt] = meanwt(x,w)
%
%         x =       vector of values to be averaged.
%         w =       corresponding weights.
%         ------------------------------------------
%         wt_mean = weighted mean.
%         wt_var =  weighted variance.
%         wt_skew = weighted skewness.
%         wt_kurt = weighted kurtosis.
%

% RE Strauss, 9/23/96
%   12/14/00 - convert input matrix to column vectors.
%   7/18/01 -  added skewness and kurtosis measures.

function [wt_mean,wt_var,wt_skew,wt_kurt] = meanwt(x,w)
  x = x(:);
  w = w(:);

  if (length(x) ~= length(w))
    error('  MEANWT: Vectors not of same length');
  end;

  w = w / sum(w);                       % Convert weights to relative freqs

  wt_mean = sum(w.*x);
  d = (x - wt_mean);

  wt_var = sum(w.*d.^2);
  wt_skew = sum(w.*d.^3);
  wt_kurt = sum(w.*d.^4);

  return;

