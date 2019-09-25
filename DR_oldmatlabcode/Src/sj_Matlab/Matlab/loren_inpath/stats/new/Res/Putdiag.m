% PUTDIAG:  Puts a vector on the diagonal of a square matrix, or matrix having 
%           more columns than rows.
%
%     Usage: matout = putdiag(matin,diag)
%
%         matin = [n x m] matrix (n <= m).
%         diag =  vector (length n), or scalar to be propagated.
%         ------------------------------------------------------------------
%         matout = corresponding [n x m] matrix with 'diag' on the diagonal.
%

% RE Strauss, 12/21/99

function matout = putdiag(matin,diag)
  [n,m] = size(matin);
  if (n>m)
    error('  PUTDIAG: input matrix is wrong size');
  end;
  
  [r,c] = size(diag);
  if (min([r,c])>1)
    error('  PUTDIAG: "diag" must be a vector or scalar');    
  end;

  if (max([r,c])==1)                    % Convert scalar to vector
    diag = diag*ones(n,1);
  end;
  
  matout = matin;
  for i = 1:n
    matout(i,i) = diag(i);
  end;

  return;
