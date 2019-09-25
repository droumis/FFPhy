% IsMonotonic: Determines whether the columns of a matrix are monotonically 
%              increasing or decreasing.
%
%     Usage: [ismono,dir] = ismonotonic(X,{alloweq})
%
%         X =       [r x c] matrix; a single row vector is transposed to a column
%                     vector.
%         alloweq = optional boolean flag indicating, if true, that consecutive 
%                     elements can be equal (i.e., not strictly monotonically 
%                     increasing or decreasing) [default = 0].
%         -----------------------------------------------------------------------
%         ismono =  [1 x c] boolean vector indicating whether columns are 
%                     monotonic.
%         dir =     [1 x c] vector containing, for each column:
%                     -1 = monotonically decreasing
%                      0 = non-monotonic
%                     +1 = monotonically increasing
%

% RE Strauss, 5/4/01

function [ismono,dir] = ismonotonic(X,alloweq)
  if (nargin < 2)  alloweq = []; end;

  if (isempty(alloweq))
    alloweq = 0;
  end;

  [isvect,ncells,iscol] = isvector(X);
  if (isvect & ~iscol)
    X = X';
  end;

  [r,c] = size(X);
  ismono = zeros(1,c);
  dir = zeros(1,c);

  for ic = 1:c
    d = X(2:r,ic) - X(1:r-1,ic);
    if (alloweq)
      if (all(d>=0))
        ismono = 1;
        dir = 1;
      elseif (all(d<=0))
        ismono = 1;
        dir = -1;
      end;
    else
      if (all(d>0))
        ismono = 1;
        dir = 1;
      elseif (all(d<0))
        ismono = 1;
        dir = -1;
      end;
    end;
  end;

  return;
