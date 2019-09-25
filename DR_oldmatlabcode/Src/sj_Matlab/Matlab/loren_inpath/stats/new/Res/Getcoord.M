% GETCOORD: Given a set of digitized point coordinates along with low & high 
%           values on the X and Y axes, recovers the original point 
%           coordinates from a scatterplot.
%
%       Usage: crds = getcoord(digcrds,axis_values,{axis_types})
%
%           digcrds =     [N+4 x 2] matrix of digitized point coordinates:
%                           1: low X tick mark
%                           2: high X tick mark
%                           3: low Y tick mark
%                           4: high Y tick mark
%                           5 to N+4: coordinates of scatter points.
%           axis_values = vector of [xlow,xhigh,ylow,yhigh].
%           axis_type =   optional vector [xtype ytype]: 0 for linear,
%                           1 for logarithmic (default [0,0]).
%           --------------------------------------------------------------
%           crds =        [N x 2] point coordinates in  scatterplot units.
%

% RE Strauss, 8/12/95
%   9/27/99 - changes in handling of input arguments.
%   2/25/00 - major rewrite; previous version used image rather than coordinate 
%               input.

function crds = getcoord(digcrds,axis_values,axis_types)
  if (nargin < 2) axis_values = []; end;
  if (nargin < 3) axis_types = []; end;

  if (isempty(axis_values))
    error('  GETCOORD: required input arguments not passed');
  end;

  if (isempty(axis_types))
    axis_types = [0 0];
  end;

  xlow = axis_values(1);
  xhigh = axis_values(2);
  ylow = axis_values(3);
  yhigh = axis_values(4);

  digcrds = rotate(digcrds,pi/4);     % Rotate entire configuration by 45 degrees
  digcrds(:,1) = digcrds(:,1) + min(digcrds(:,1));
  digcrds(:,2) = digcrds(:,2) + min(digcrds(:,2));

  ticks = digcrds(1:4,:);
  x = digcrds(5:end,1);
  y = digcrds(5:end,2);

  xlog = 0;
  ylog = 0;
  if (axis_types(1))
    xlog = 1;
    xlow = log(xlow);
    xhigh = log(xhigh);
  end;
  if (axis_types(2))
    ylog = 1;
    ylow = log(ylow);
    yhigh = log(yhigh);
  end;

  % Get linear parameters of X,Y axes from low and high points

  b = lineeqn(ticks(1,:),ticks(2,:));
  xintcp = b(1);
  xslope = b(2);
  xscale = (xhigh-xlow)/eucl(ticks(1:2,:));

  b = lineeqn(ticks(3,:),ticks(4,:));
  yintcp = b(1);
  yslope = b(2);
  yscale = (yhigh-ylow)/eucl(ticks(3:4,:));

  % Convert other points to cartesian coordinates of graph

  xscores = project(digcrds(:,1),digcrds(:,2),xslope,xintcp);
  yscores = project(digcrds(:,1),digcrds(:,2),yslope,yintcp);

  xscores = (xscores - xscores(1)) ./ (xscores(2)-xscores(1));
  yscores = (yscores - yscores(3)) ./ (yscores(4)-yscores(3));

  xscores = xscores*(xhigh-xlow) + xlow;
  yscores = yscores*(yhigh-ylow) + ylow;

  if (xlog)
    xscores = exp(xscores);
  end;
  if (ylog)
    yscores = exp(yscores);
  end;

  crds = [xscores(5:end) yscores(5:end)];

  return;
