% CORRCMPL:  Given a data matrix, returns a symmetric distance matrix, where 
%            the distances are the complements of the rank correlations among
%            columns.
%
%     Usage: dist = corrcmpl(X)
%
%            X = (n x p) data matrix.
%            ------------------------------------------------
%            dist = (n x n) square symmetric distance matrix.
%

% RE Strauss, 11/6/97
%   9/16/99 - put exact zeros on diagonal.

function dist = corrcmpl(X)
  dist = rankcorr(X);
  dist = 1-dist;
  dist = putdiag(dist,0);

  return;

