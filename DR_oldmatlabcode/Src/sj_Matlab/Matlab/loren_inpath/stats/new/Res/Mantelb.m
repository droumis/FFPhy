% MANTELB:  Returns a bootstrapped sample from a distance matrix based on a 
%           parametric bootstrap, given the standard errors of the distances.
%
%     Usage: Db = mantelb(D,Dstderr)
%
%           D =       square symmetric distance matrix.
%           Dstderr = corresponding matrix of standard errors.
%           --------------------------------------------------
%           Db =      bootstrapped distance matrix.
%

% RE Strauss, 10/19/00

function Db = mantelb(D,Dstderr)
  d = trilow(D);                          % Convert to col vectors
  ds = trilow(Dstderr);

  Db = trisqmat(d + randn(size(d)).*ds);  % Add noise to observed distances

  return;
