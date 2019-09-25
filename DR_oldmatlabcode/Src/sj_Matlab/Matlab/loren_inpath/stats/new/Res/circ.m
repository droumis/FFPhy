% CIRC: Return points on the circumference of a circle, and plot the circle.
%
%     Usage: [x,y] = circ(radius,{center},{noplot},{npts})
%
%       radius = scalar value for radius of circle.
%       center = optional vector (length 2) of x,y coordinates of center of circle 
%                   [default = (0,0)].
%       noplot = optional boolean value indicating, if true, that the plot is to 
%                   be suppressed [default = 0].
%       npts =   optional number of points returned [default = 100].
%

function [x,y] = circ(radius,center,noplot,npts)

  if (nargin < 1)
    error('CIRC: radius not passed');
  end;

  if (nargin < 2)
    center = [];
  end;
  if (nargin < 3)
    noplot = [];
  end;
  if (nargin < 4)
    npts = [];
  end;

  if (isempty(center))
    center = [0 0];
  end;
  if (isempty(noplot))
    noplot = 0;
  end;
  if (isempty(npts))
    npts = 100;
  end;

  theta = linspace(0,2*pi,npts);

  x = radius .* cos(theta) + center(1);
  y = radius .* sin(theta) + center(2);

  if (~noplot)
    plot(x,y);
     hold on;
    plot(x,y,'x');
    plot(center(1),center(2),'d');
    hold off;
    axis square;
  end;

  return;