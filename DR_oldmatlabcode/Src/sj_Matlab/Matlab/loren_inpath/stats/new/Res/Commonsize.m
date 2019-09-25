% COMMONSIZE: Checks whether two matrices are identical in size.  If not, and if 
%             one is a scalar, propagates its value to the size of the other.
%
%     Usage: [samesize,x,y] = commonsize(x,y,{c})
%
%         x,y =       matrices of identical size, or matrix and a scalar.
%         c =         optional matrix for size adjustment, if both x and y are 
%                       scalars [default = null].
%         --------------------------------------------------------------------
%         samesize =  boolean vector indicating, if true, that output matrices 
%                       are compatible.
%         x,y =       input matrices, expanded if a scalar.
%

% RE Strauss, 6/7/00

function [samesize,x,y] = commonsize(x,y,c)
  if (nargin < 3) c = []; end;

  samesize = 1;
  sizex = size(x);
  sizey = size(y);
  isscalx = isscalar(x);
  isscaly = isscalar(y);

  if (any(sizex~=sizey) | (isscalx & isscaly))

    if (isscalx & ~isscaly)
      x = x * ones(size(y));
    elseif (~isscalx & isscaly)
      y = y * ones(size(x));
    elseif (isscalx & isscaly)
      x = x * ones(size(c));
      y = y * ones(size(c));
      if (isempty(c))
        samesize = 0;
      end;
    else
      samesize = 0;
    end;
  end;

  return;
