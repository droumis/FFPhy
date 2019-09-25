% LOGISTBT: Objective LOGISTIC function for bootstrapped confidence intervals.
%
%     Usage: p = func([t,W])
%
%           t = vector of abscissa (time) values.
%           W = vector of ordinate (size) values.
%           -------------------------------------
%           p = current parameter estimates.
%

% RE Strauss, 4/8/99
%   1/4/00 -    changed fminu() to fmins().

function p = func(X,notused1,notused2,notused3)
  t = X(:,1);
  W = X(:,2);

  k = 1;                                      % Initial parameter estimates
  h = mean(t);                    
  Wa = max(W);

  p = fmins('logistf',[k h Wa],[],[],t,W);    % Fit model

  return;
