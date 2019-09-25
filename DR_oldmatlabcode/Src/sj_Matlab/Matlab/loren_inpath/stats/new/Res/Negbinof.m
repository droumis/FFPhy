% NEGBINOF: Objective function to determine the maximum-likelihood estimator of 
%           the negative-binomial parameter k.  See Young & Young (1998): 7-8.
%
%     Usage: d = negbinof(k,mu,P)
%
%           k =   current estimate of k.
%           mu =  estimate of parameter mu.
%           P =   empirical pdf for X = 0:maxcount.
%           -------------------------------------------------
%           d =   inequality difference for the MLE equation.
%
%   <<< Not working: allows k to become arbitrarily large. >>> 
%

% RE Strauss, 7/1/99

function d = negbinof(k,mu,P)
  L = sum(P) * log(1 + mu/k);         % Left term in MLE equation

  R = 0;                              % Right term
  t = 0;
  for i = 2:length(P)
    j = i-1;
    t = t + 1./(k+j-1);
    R = R + P(i)*t;
  end;

  d = abs(L-R);                       % Magnitude of inequality
%k_L_R_d = [k L R d]

  return;