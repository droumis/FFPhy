% BINMATNO: Estimates number of possible (0,1)-element [binary] matrices given
%           the row and column sums, based on Eqn 3.3 of Hendrickson (1995).
%           The estimate is usually within 10% of the true value and,
%           for matrices with many rows, often accurate to within 2%.
%
%     Syntax: A = binmatno(rowsums,colsums)
%
%           rowsums = vector of row sums.
%           colsums = vector of column sums.
%           -------------------------------------------------------
%           A = estimated number of possible binary matrices.
%

function A = binmatno(rowsums,colsums)
  N = sum(rowsums);
  if (N ~= sum(colsums))
    error('  Totals of row-sum and column-sum vectors must be identical (N)');
  end;

  r = length(rowsums);
  s = length(colsums);
  f = zeros(1,s+1);
  A = 0;

  for j = 0:s
    f(j+1) = sum(rowsums==j);
  end;

  j = 0:s;
  sumfjsj = sum(f.*j.*(s-j));
  if (sumfjsj==0)
    return;
  end;

  term1 = 2*sumfjsj/s;

  term2 = ((s-1)/(pi.*term1))^((s-1)/2);

  term3 = 1;
  for j = 0:s
    sfj = combin(s,j)^f(j+1);
    term3 = term3*sfj;
  end;

  term4 = exp(-(s-1)*(sum(colsums.^2)-(N*N/s))/term1);

  A = term3*term2*sqrt(s)*term4;

  return;
