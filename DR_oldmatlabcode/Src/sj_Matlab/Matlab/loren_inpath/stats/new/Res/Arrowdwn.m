% ARROWDWN: Plots downard-pointing arrow(s) to mark an x-axis position on 
%           a histogram or other plot.
%
%     Usage: arrowdwn(x,y,alength)
%
%           x,y =     absolute coordinates of bottom (tip) of arrow.  If x,y 
%                       are vectors of same length, then one arrow is plotted 
%                       for each pair of coordinates.
%           alength = arrow length in y-units [default = 0.15*length of y-axis].
%

% RE Strauss, 3/27/97
%   8/19/99 - make arrow length absolute rather than relative to y-axis length.

function arrowdwn(x,y,alength)
  if (nargin < 3) alength = []; end;

  if (length(x) ~= length(y))
    error('  ARROWDWN: Coordinate vectors must be same length');
  end;

  v = axis;
  xrng = v(2)-v(1);
  yrng = v(4)-v(3);

  if (isempty(alength))
    alength = 0.15 * yrng;
  end;
  head = 0.2;

  shaft = alength * yrng;
  hy = head * shaft;
  hx = 0.4 * head * alength * xrng;

  hold on;
  for a = 1:length(x)
    x1 = x(a);
    y1 = y(a);
    x2 = x1;
    y2 = y1 + shaft;
    x3 = x1 - hx;
    y3 = y1 + hy;
    x4 = x1 + hx;
    y4 = y3;
    plot([x2 x1 x3 x1 x4],[y2 y1 y3 y1 y4],'k');
  end;

  return;



