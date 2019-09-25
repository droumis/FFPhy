% NSCORES:  Estimates normal scores for a sample of size n.  See also rankits().
%
%     Usage: scores = nscores(n)
%

% RE Strauss, 12/8/02 - completely rewritten from previous version.

function scores = nscores(n)

  if (~isintegr(n) | n<=0)
    error('  NSCORES: n must be a positive integer.');
  end;

  scores = norminv(1/(n+1):1/(n+1):n/(n+1))';     % Probability scores
  
  return;
  

    