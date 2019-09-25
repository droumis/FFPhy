% SUITES: Evaluates existence in a data matrix of character (variable) suites, 
%         subsets of correlated characters.
%
%     Usage: suites(X)
%
%         X = [n x p] data matrix.
%

% RE Strauss, 9/17/99

function suites(X)
  [n,p] = size(X);

  C = corrcoef(X);                      % Correlation matrix
  dist = 1-abs(C);                      % Distance matrix among variables

  topology = upgma(dist);
  v = axis;
  v(1:2) = [0 1];
  axis(v);
  puttick(0:0.1:1);

  return;
