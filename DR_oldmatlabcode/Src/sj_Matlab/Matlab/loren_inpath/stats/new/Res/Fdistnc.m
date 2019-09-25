% FDISTNC: Returns cumulative probability of f for the noncentral
%          F distribution with parameters v1,v2 and noncentrality lambda.
%
%     Syntax: prob = fdistnc(f,v1,v2,lambda)
%

% Lenth, R.V. 1987.  Algorithm AS 226: Computing noncentral beta probabilities.
%   Appl. Stat. 36:241-244.
% Frick, H. 1990. Algorithm AS R84: A remark on Algorithm AS 226: Computing
%   noncentral beta probabilities.  Appl. Stat. 39:311-312.

function prob = fdistnc(f,v1,v2,lambda)
  prob = NaN;

  x = (v1.*f)./(v1.*f+v2);
  prob = betanc(x,v1/2,v2/2,lambda);

  return;

