% PLOTCIRCLE: plots a circle
%
%     Usage: plotcircle({radius},{n})
%
%     radius = optional radius of circle [default = 1].
%     n =      optional number of plotted points [default = 100].
%

function plotcircle(radius,n)
  if (nargin < 1)
    radius = [];
  end;
  if (nargin < 2)
    n = [];
  end;
  
  default_radius = 1;                 % Default input argument values
  default_n = 100;
  
  if (isempty(radius))
    radius = default_radius;
  end;
  if (isempty(n))
    n = default_n;
  end;

  theta = linspace(0,2*pi,n);
  
  x = radius .* cos(theta);
  y = radius .* sin(theta);
  
  plot(x,y,'k');
  axis square;
  
  return;
  