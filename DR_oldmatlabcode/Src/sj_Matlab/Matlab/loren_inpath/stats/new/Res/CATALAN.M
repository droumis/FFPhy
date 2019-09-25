% CATALAN: Calculates a column vector of the nth Catalan numbers, where
%         n varies from n1 to n2 and n2>n1.
%         If only a single integer argument is supplied, the corresponding
%         Catalan number is returned as a scalar.
%
%     Usage: cat = catalan(n1,n2)
%

function cat = catalan(n1,n2)
  if (nargin<2)
    cat = (1/(n1+1))*combin(2*n1,n1);
  else
    if (n1>n2)
      error('  Error: input arguments must be in increasing sequence');
    end;

    nval = n2-n1+1;
    cat = zeros(nval,1);
    for i=1:nval
      j = i+n1-1;
      cat(i) = (1/(j+1))*combin(2*j,j);
    end;
  end;

  return;
