% ZSCORE: Standardizes the columns of an [n,p] matrix X.  Ignores missing data.
%
%     Usage: z = zscore(X)
%
%         X = [n x p] data matrix.
%         ----------------------------------------
%         z = [n x p] matrix of standardized data.
%

% RE Strauss, 10/14/97
%   6/22/00 - allow for columns of constants; allow for missing data.

function z = zscore(X)
  if (nargin < 1) help zscore; return; end;

  stdX = std(X);
  meanX = mean(X);
  ismiss = misscheck(stdX);
  if (ismiss)
    stats = univar(X,[1 0 0 1]);
    meanX = stats(1,:);
    stdX = stats(2,:);
  end;

  const_col = 0;
  if (any(stdX<eps))
    i = find(stdX<eps);
    stdX(i) = NaN*ones(1,length(i));
    const_col = 1;
  end;

  [n,p] = size(X);
  e = ones(n,1);
  X = X - e*meanX;
  z = X ./ (e*stdX);

  if (const_col)
    z(:,i) = zeros(n,length(i));
  end;
  if (ismiss)
    [i,j] = find(~isfinite(X));
    for k = 1:length(i)
      z(i(k),j(k)) = NaN;
    end;
  end;

  return;
