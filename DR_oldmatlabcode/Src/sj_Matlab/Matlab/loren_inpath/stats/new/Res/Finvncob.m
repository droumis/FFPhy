% FINVNCOB: Objective function for FINVNC
%
%     Syntax: delta = finvncob(x,prob,v1,v2,lambda)
%
%         x =      current value of noncentral beta distribution
%         prob =   target cumulative probability
%         v1,v2 =  degrees of freedom
%         lambda = noncentrality parameter
%         ======================================================
%         delta =  abs(cdf(F) - prob)
%

% RE Strauss, 4/7/96

function delta = finvncob(f,prob,v1,v2,lambda)
  p = fdistnc(f,v1,v2,lambda);
  delta = abs(p-prob);

  return;
