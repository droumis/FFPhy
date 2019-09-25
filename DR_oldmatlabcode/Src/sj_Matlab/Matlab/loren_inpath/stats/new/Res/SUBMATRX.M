% SUBMATRX: Part of a Wright-style factor analysis.
%           Given a square correlation matrix and a vector of indices to 
%           variables (not necessarily contiguous), extracts and returns the 
%           submatrix.
%
%     Usage: sub = submatrx(matrx,vars)
%
%           matrx = [n x n] correlation matrix
%           vars =  m-length vector of vars (specified by row=column)
%           sub =   [m x m] submatrix
%

% RE Strauss, 1/21/95

function sub = submatrx(matrx,vars)
  sub = matrx(vars,vars);
  return;

