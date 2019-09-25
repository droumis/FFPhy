% ZScores: Standardizes the columns of an [n,p] matrix X.  Ignores missing data.
%
%     Usage: z = zscores(X)
%
%         X = [n x p] data matrix.
%         ----------------------------------------
%         z = [n x p] matrix of standardized data.
%

% RE Strauss, 1/3/03

function z = zscores(X)
  z = zscore(X);
  return;
  