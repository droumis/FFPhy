% TREEVECT: Given the number of terminal taxa T and a value of Rohlf's (1983) 
%           integer M, recovers the N-tuple N characterizing the tree topology.
%
%       Usage: N = treevect(T,M)
%
%           T = number of terminal taxa (assumed to be in fixed sequence).
%           M = Rohlf's topology number, from 0 to treenum(T)-1.
%           --------------------------------------------------------------
%           N = row vector characterizing tree topology.
%

% Rohlf, FJ. 1983. Numbering binary trees with labeled terminal vertices.  Bull. Math. 
%   Biol. 45:433-40.

function N = treevect(T,M)
  ntree = treenum(T);
  if ((M < 0) | (M > ntree-1))
    disp(sprintf('  Error: Allowable range of M = 0-%1.0f', ntree-1));
    error(' ');
  end;

  N = zeros(1,T);

  A = M;
  for i = 3:T-1
    B = floor(A/(2*i-3));
    N(i) = A - (2*i-3)*B;
    A = B;
  end;
  N(T) = A;
  return;
