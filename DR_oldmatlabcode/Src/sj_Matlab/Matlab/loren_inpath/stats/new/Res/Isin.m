% ISIN: Returns a boolean matrix for matrix A, indicating whether 
%       or not the entries of A are contained within the matrix B.  
%       Matrices A and B do not have to be the same size, but can have a maximum 
%       of two dimensions.
%
%     Usage: [b,i,j] = isin(A,B,{tol})
%
%       A =   test matrix.
%       B =   target matrix.
%       tol = optional equality tolerance [default = eps].
%       ---------------------------------------------------------------
%       b =   boolean matrix corresponding in size to matrix A.
%       i =   row index of B containing first occurrence of value in A.
%       j =   col index of B containing first occurrence of value in A.
%

% RE Strauss, 11/10/98
%   5/31/01 - add option of output indices.

function [b,ind1,ind2] = isin(A,B,tol)
  if (nargin < 3) tol = []; end;

  get_index = 0;
  if (nargout > 1)
    get_index = 1;
  end;

  if (isempty(tol))
    tol = eps;
  end;
  
  [r,c] = size(A);
  b = zeros(size(A));
  if (get_index);
    ind1 = zeros(size(A));
    ind2 = zeros(size(A));
  end;
  Bv = B(:);

  for i = 1:r
    for j = 1:c
      if sum(abs(A(i,j)-Bv)<(tol+eps))
        b(i,j) = 1;
        if (get_index)
          bb = abs(A(i,j)-B)<(tol+eps);
          [ii,jj] = find(bb);
          ind1(i,j) = ii(1);
          ind2(i,j) = jj(1);
        end;
      end;
    end;
  end;

  if (get_index)
    if (~ismatrix(B))
      if (sum(ind1)<sum(ind2))
        ind1 = ind2;
      end;
      ind2 = [];
    end;
  end;

  return;
