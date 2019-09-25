% Iscorr: Returns 1 if the input matrix is in the form of a correlation matrix, 
%         and 0 if not.  Doesn't check for valid combinations of correlations,
%         or for missing data.
%
%     Usage: b = iscorr(c)
%

% RE Strauss, 2/8/00

function b = iscorr(c)
  [r,p] = size(c);
  b = 1;

  if (r~=p)
    b = 0;
    return;
  end;

  if (sum(diag(c))~=r)
    b = 0;
  end;
  if (sum(trilow(c)-trilow(c')) > eps*r*(r-1))
    b = 0;
  end;

  return;
