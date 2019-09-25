% BIMODALF: Objective function for BIMODAL.  Coefficient is centered on zero,
%           the value for a uniform distribution.  Positive values are
%           increasingly bimodal, negative values are increasingly unimodal.
%

% RE Strauss, 10/21/97

function b = bimodalf(X,grps,initsoln,nulldist,nullflag)
  [n,p] = size(X);

  if (nargin < 4)
    nulldist = 0;
  end;

  if (nulldist)                     % Random-uniform sample
    X = rand(n,1)*ones(1,p);          % Identical sample for all variables
  end;

  kurt = kurtosis(X);               % Kurtosis, by column
  b = log(1.8./kurt);               % Bimodality coeff, by column
%  bias = 0.0294-(0.4051./sqrt(n));  % Bias correction
%  b = b - bias;

  return;

