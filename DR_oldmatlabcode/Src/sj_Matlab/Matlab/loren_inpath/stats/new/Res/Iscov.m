% Iscov: Returns 1 if the input matrix is in the form of a covariance matrix, 
%        and 0 if not.  Doesn't check for valid combinations of covariances,
%        or for missing data.
%
%     Usage: b = iscov(c)
%

% RE Strauss, 6/7/00 (modified from iscorr.m)

function b = iscov(c)
  [r,p] = size(c);
  b = 1;

  if (r~=p)
    b = 0;
    return;
  end;

  if (sum(trilow(c)-trilow(c')) > eps*r*(r-1))
    b = 0;
  end;

  return;
