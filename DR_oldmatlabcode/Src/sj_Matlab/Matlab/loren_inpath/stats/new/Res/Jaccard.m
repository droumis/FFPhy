% JACCARD: Calculates the Jaccard measure of co-occurrence.
%
%     Syntax:  J = jaccard(distrib)
%
%           distrib = [S x N] binary occurrence matrix for N species 
%                       across S sites.
%           ---------------------------------------------------------------------
%           J =       [N x N] symmetric matrix of coefficient values.
%
%     See jaccardd() for the corresponding distance measure.
%

% RE Strauss, 2/19/97

function J = jaccard(distrib)
  i = (distrib(:)~=0 & distrib(:)~=1);
  if (sum(i)>0)
    error('  JACCARD: Binary (presence/absence) data only');
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

  return;
