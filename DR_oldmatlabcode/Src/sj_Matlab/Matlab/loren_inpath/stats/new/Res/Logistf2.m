% LOGISTF2: Objective function for LOGISTIC for 2-parameter model.
%
%     Syntax: ypred = logistf2(b,x)
%
%           b = current parameter estimates.
%           x = vector of abscissa values.
%           -------------------------------------
%           ypred = predicted ordinate values.
%

function ypred = logistf2(b,x)
  ypred = 1./(1+exp(b(1)-b(2).*x));
  return;

