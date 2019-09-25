% LOGISTF3: Objective function for LOGISTIC for 3-parameter model.
%
%     Syntax: ypred = logistf3(b,x)
%
%           b = current parameter estimates.
%           x = vector of abscissa values.
%           -------------------------------------
%           ypred = predicted ordinate values.
%

function ypred = logistf3(b,x)
  ypred = b(3)./(1+exp(b(1)-b(2).*x));
  return;

