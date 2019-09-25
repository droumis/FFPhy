% IsOutOfBounds: determines whether subscripts are out of bounds for a particular matrix.
%
%     Usage: out = isoutofbounds(X,s)
%
%         X = [n1 x n2 x ... x np] P-dimensional matrix.
%         s = [m x p] matrix of subscripts, where each row indexes a cell in X.
%         ---------------------------------------------------------------------
%         out = [m x p] boolean vector indicating subscripts out of bounds.
%

% RE Strauss, 5/14/03

function out = isoutofbounds(X,s)
  sizeX = size(X);
  dim = ndims(X);
  
  [m,p] = size(s);
  if (p~=dim)
    error('  IsOutOfBounds: number of cols of subscript matrix must match dimensions of X.');
  end;
  
  out = zeros(m,p);
  
  for im = 1:m
    i = find(s(im,:)<1 | s(im,:)>sizeX);
    if (~isempty(i))
      out(im,i) = ones(1,length(i));
    end;
  end;

  return;
  