% SUMSQDIFFF: Objective function for sumsqdiff()
%
%     Usage: dist = sumsqdiff(X,nu1,nu2,nu3)
%
%         X =     [r x c] matrix, for which sum-of-squared differences are 
%                   calculated among columns, assuming that rows have been 
%                   bootstrapped.
%         nu.. =  arguments passed by bootstrp, but unused.
%         ----------------------------------------------------------------
%         dist =  [1 x c] vector of pairwise distances among cols of X; 
%                   concatenated rows of upper diagonal.
%

% RE Strauss, 5/16/00

function dist = sumsqdiff(X,nu1,nu2,nu3)
  [r,c] = size(X);
  dist = zeros(c,c);
  
  for i = 1:(c-1)                         % Distance matrix
    for j = (i+1):c
      d = X(:,i) - X(:,j);
      dist(j,i) = d'*d;
    end;
  end;

  dist = trilow(dist)';

  return;
