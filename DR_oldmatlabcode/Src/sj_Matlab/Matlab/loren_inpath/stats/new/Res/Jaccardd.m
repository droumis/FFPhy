% JACCARDD: Calculates the Jaccard distance measure.
%
%     Syntax:  dist = jaccardd(distrib)
%
%           distrib = [S x N] binary occurrence matrix for N species 
%                       across S sites.
%           ---------------------------------------------------------------------
%           J =       [N x N] symmetric distance matrix.
%
%     See jaccard() for the corresponding measure of association.
%

function dist = jaccardd(distrib)
  i = (distrib(:)~=0 & distrib(:)~=1);
  if (sum(i)>0)
    error('  JACCARDD:  Binary (presence/abs ence) data only');
  end;

  [S,N] = size(distrib);
  J = ones(N,N);

  for i = 1:(N-1)
    for j = (i+1):N
      joint = sum(distrib(:,i) & distrib(:,j));
      either = sum(distrib(:,i) | distrib(:,j));
      if (either > 0)
        J(i,j) = joint / either;
      else
        J(i,j) = 0;
      end;
      J(j,i) = J(i,j);
    end;
  end;

  dist = 1-J;

  return;
