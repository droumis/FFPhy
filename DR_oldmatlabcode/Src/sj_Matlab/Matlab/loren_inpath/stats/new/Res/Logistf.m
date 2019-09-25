 % LOGISTF: Objective function for LOGISTIC.
%
%     Syntax: mse = logistf(p,t,W)
%
%           p = current parameter estimates.
%           t = vector of abscissa (time) values.
%           W = vector of ordinate (size) values.
%           -------------------------------------
%           mse = mean squared error.
%

% RE Strauss, 4/7/99

function mse = logistf(p,t,W)
  k = p(1);
  h = p(2);
  Wa = p(3);

  Wpred = Wa./(1+exp(-k*(t-h)));
  mse = (W-Wpred)'*(W-Wpred);

  return;

