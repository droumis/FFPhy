% RandPermCols: Randomly permute columns of a matrix
%
%     Usage: Y = randpermcols(X)
%

% RE Strauss, 2/7/02

function Y = randpermcols(X)
  [r,c] = size(X);
  Y = X;
  
  for i = 1:c
    Y(:,i) = Y(randperm(r),i);
  end;

  return;