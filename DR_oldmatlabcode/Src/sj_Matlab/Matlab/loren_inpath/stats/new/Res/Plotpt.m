% PLOTPT: Add a single point to an existing plot, extending axis ranges if 
%         necessary.
%
%     Usage: plotpt(crds,{symbol})
%
%         crds =    2-element vector of cartesian coordinates.
%         symbol =  character string indicating symbol/color for plotting 
%                     [default = 'ko'].
%

% RE Strauss, 1/4/00

function plotpt(crds,symbol)
  if (nargin < 2) symbol = []; end;

  if (isempty(symbol))
    symbol = 'ko';
  end;

  if (length(crds)~=2)
    error('  PLOTPT: input coordinates must be given as 2-element vector');
  end;

  x = crds(1);
  y = crds(2);

  hold on;
    plot(x,y,symbol);
  hold off;

  v = axis;                             % Get current axis bounds
  vs = v;
  d = 0.02;                             % Edge buffer for new point

  if (x < v(1))                         % If point out of range, extend bounds
    v(1) = x - d*(v(2)-v(1));
  elseif (x > v(2))
    v(2) = x + d*(v(2)-v(1));
  end;
  if (y < v(3))
    v(3) = y - d*(v(4)-v(3));
  elseif (y > v(4))
    v(4) = y + d*(v(4)-v(3));
  end;
  if (any(v ~= vs))
    axis(v);
  end;

  return;