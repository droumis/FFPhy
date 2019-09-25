% Partf:  Recursive function for partitioning N linearly spaced objects into k 
%         groups.  See partion().
%
%     Usage: grps = partf(N,k,grps)
%
%             N = number of objects.
%             k = number of groups.
%             ----------------------------------------
%             grps = [npart x N] matrix of partitions.
%

% RE Strauss, 8/27/99

function grps = partf(N,k,grps)
  if (nargin < 3) grps = []; end;

  if (k==1)
    grps = ones(1,N);
    return;
  end;

  for i = 1:(N-k+1)
    g = partf(N-i,k-1);
    [r,c] = size(g);
    grps = [grps; k*ones(r,i) g];
  end;

  return;
