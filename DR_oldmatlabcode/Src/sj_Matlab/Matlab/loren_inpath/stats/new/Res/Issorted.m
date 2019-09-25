% ISSORTED: Determines whether the columns of a matrix are in sequence.
%
%     Usage: b = issorted(X,{dir})
%
%           X =   [r x c] input matrix.
%           dir = optional value indicating whether the sequence should be 
%                   increasing (zero or positive value, the default) or 
%                   decreasing (negative value).
%           ------------------------------------------------------------------
%           b =   [1 x c] boolean vector, with entries of 1 for columns that 
%                 are in sequence, and 0 for columns that are not in sequence.
%

% RE Strauss, 2/19/00

function b = issorted(X,dir)
  if (nargin < 2) dir = []; end;

  if (isempty(X))
    b = [];
    return;
  end;
  if (isempty(dir))
    dir = 0;
  end;

  if (dir(1)<0)
    X = -X;
  end;

  [r,c] = size(X);
  b = ones(1,c);

  if (r>1)
    for j = 1:c
      if (any(X(2:r,j) < X(1:r-1,j)))
        b(j) = 0;
      end;
    end;
  end;

  return;
