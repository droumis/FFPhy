% CORRADJ:  Creates a first-order adjacency matrix for a square matrix,
%           containing 1's for neighboring cell combinations and 0's
%           elsewhere.  The adjacency matrix is of order N*N and indexes 
%           the concatenated rows (or columns) of an associated square matrix.
%
%       Usage: A = corradj(n)
%
%           n = order of corresponding square matrix matrix
%           A = sparse [n^2 x n^2] adjacency matrix
%

% RE Strauss, 2/15/95

function w = corradj(n)
  n2 = n*n;                       % Number of matrix elements

  % Construct adjacency matrix
  w = zeros(n2,n2);
  i = 0;
  for r1=1:n
    for c1=1:n
      i = i+1;
      for r2=r1:min(r1+1,n)
        if (all([r2==r1,c1<n]))
          j = i+1;
          w(i,j)=1;
          w(j,i)=1;
        elseif (r2>r1)
          for c2=max(c1-1,1):min(c1+1,n)
            j = (r2-1)*n+c2;
            w(i,j)=1;
            w(j,i)=1;
          end;
        end;
      end;
    end;
  end;
  w = sparse(w);
  return;


