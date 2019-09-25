% SUMSQSCALE: Scales the columns of a matrix so that the squared elements sum to 
%             unity.
%
%     Usage: S = sumsqscale(X)
%
%             X = [r,c] matrix.
%             --------------------------------
%             S = corresponding scaled matrix.
%

% RE Strauss, 5/2/00

function S = sumsqscale(X)
  [r,c] = size(X);

  s = sign(X);
  sq = X.^2;
  ssq = sum(sq);
  for i = 1:c
    sq(:,i) = sq(:,i)./ssq(i);
  end;
  S = s.*sqrt(sq);

  return;
