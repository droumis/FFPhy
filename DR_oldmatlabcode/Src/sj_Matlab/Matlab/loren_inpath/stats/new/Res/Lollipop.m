% LOLLIPOP: Produces 3D scatterplot in which points are anchored onto the 
%           X,Y plane by a stem (subtended line).  Does not open a new figure window.
%
%     Syntax:  lollipop(x,y,z,{color})  OR  lollipop([x,y,z],{color})
%
%           x,y,z = column vectors of identical length.
%           color = optional character value indicating color or circle and line
%                     [default = 'k'].
%

% RE Strauss, 2/18/98
%   11/25/99 - changed plot colors for Matlab v5.
%   8/23/01 -  added box statement for Matlab v6.
%   6/6/03 -   added optional color indicator; use putbnd3() to set axis bounds.

function lollipop(x,y,z,color)
  if (nargin < 3)
    if (nargin < 2)
      color = [];
    else
      color = y;
    end;
    z = x(:,3);
    y = x(:,2);      
    x = x(:,1);
  elseif (nargin < 4)
    color = [];
  end;
  
  if (isempty(color)) color = 'k'; end;

  if (length(x)~=length(y) | length(x)~=length(z))
    error('  LOLLIPOP: input vectors must be equal in length.');
  end;

  if (~isscalar(color))
    error('  LOLLIPOP: color must be scalar character value.');
  end;
  
  v = putbnd3([x y z],[],1);
  if (min(z)>=0 & min(z)<max(z)/2)
    v(5) = 0;
  end;
  zmin = v(5);
  
  hold on;
  plot3(x,y,z,[color,'o']);
  for i=1:length(x)
    plot3([x(i) x(i)],[y(i) y(i)],[zmin z(i)],[color,'-']);
  end;
  axis(v);
  box on;
  view(-24,28);
  hold off;

  return;

