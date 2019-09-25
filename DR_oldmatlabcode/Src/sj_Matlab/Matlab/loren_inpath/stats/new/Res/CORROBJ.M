% CORROBJ:  Calculates the objective function (Moran's autocorrelation
%           index I) for a square symmetric matrix, given the corresponding
%           adjacency matrix.
%
%     Usage: fn = corrobj(r,w)
%
%           r =  [n x n] square symmetric matrix
%           w =  [n^2 x n^2] adjacency matrix
%           fn = value of objective function
%

% RE Strauss, 3/11/95

function fn = corrobj(r,w)
  nr = size(r,1);
  nw = size(w,1);
%  if (nw~=nr*nr)
%    error('  Input and adjacency matrices incompatible');
%  end;

  sumw = sum(sum(w));
  sumz2 = sum(sum(r.*r));

  rv = r(:);
  i = ones(nw,1);
  fn = nr*(i'*(rv*rv'.*w)*i)/sumw/sumz2;
  return;


