% WRAP: Remainder after division.  Differs from the modulus function rem()
%       in that wrap(n,n) returns the value n rather than zero, and also
%       wraps a negative value around to a positive one.  Useful for
%       handling wrap-around vectors (such as polygons).
%       Returns a vector w of the same size as x.  If y is a vector
%       size-compatible with x, element-wise operations are performed.
%       If y is a scaler, it is used as the denominator for all elements of
%       x.
%
%     Usage: w = wrap(x,y)
%

% RE Strauss, 10/26/94

function w = wrap(x,y)
  nx = max(size(x));
  w = x;

  if (nx > 1)                   % If x is a vector
    if (max(size(y)) == 1)      % and y is a scaler
      y = y*ones(size(x));      % make y a vector
    end;
  end;

  for i=1:nx
    if (x(i) > 0)
      w(i) = rem(x(i),y(i));
      if (w(i) == 0)
        w(i) = y(i);
      end;
    else
      while (x(i) <= 0)
        x(i) = x(i) + y(i);
      end;
      w(i) = x(i);
    end;
  end;

  return;

