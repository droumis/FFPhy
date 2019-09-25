% GETCOORD: Given an image specified by [X,map] and low & high values on the
%           X and Y axes, permits a series of points to be digitized and
%           recovers the original point coordinates.
%
%       Usage: [x,y] = getcoord(X,map,ticks,axis_values,{axis_types})
%
%           X =           image matrix.
%           map =         colormap for image.
%           ticks =       coords of axis tick marks, from getaxes().
%           axis_values = vector of [xlow,xhigh,ylow,yhigh].
%           axis_type =   optional vector [xtype ytype]: 0 for linear,
%                           1 for logarithmic (default [0,0]).
%           ----------------------------------------------------------------
%           x,y =         point coordinates in original scatterplot units.
%

% RE Strauss, 8/12/95
%   9/27/99 - changes in handling of input arguments.

function [x,y] = getcoord(X,map,ticks,axis_values,axis_types)
  if (nargin < 5) axis_types = []; end;

  if (nargin < 4)
    error('GETCOORD: required input arguments not passed');
  end;

  if (isempty(axis_types))
    axis_types = [0 0];
  end;

  sintheta = sin(pi/4);               % For rotations of 45 degrees
  costheta = cos(pi/4);

  xlow = axis_values(1);
  xhigh = axis_values(2);
  ylow = axis_values(3);
  yhigh = axis_values(4);

  if (axis_types(1))
    xlog = 1;
    xlow = log(xlow);
    xhigh = log(xhigh);
  else
    xlog = 0;
  end;

  if (axis_types(2))
    ylog = 1;
    ylow = log(ylow);
    yhigh = log(yhigh);
  else
    ylog = 0;
  end;

  disp('  Digitize data points...');
  imshow(X,map);
  [x,y] = impixel;

  % Rotate entire configuration by 45 degrees

%  xin = x - min(x);
%  yin = -y - min(-y);
  xin = x;
  yin = -y;
  x = xin*costheta - yin*sintheta;
  y = yin*costheta + xin*sintheta;

  % Get linear parameters of X,Y axes from low and high points

  xslope = (ticks(2,2)-ticks(1,2))/(ticks(2,1)-ticks(1,1));
  xintcp = ticks(1,2) - xslope*ticks(1,1);
  xscale = (xhigh-xlow)/eucl(ticks(1:2,:));
  xlowx =  ticks(1,1);
  xlowy =  ticks(1,2);

  yslope = (ticks(4,2)-ticks(3,2))/(ticks(4,1)-ticks(3,1));
  yintcp = ticks(3,2) - yslope*ticks(3,1);
  yscale = (yhigh-ylow)/eucl(ticks(3:4,:));
  ylowx =  ticks(3,1);
  ylowy =  ticks(3,2);

  % Convert other points to cartesian coordinates of graph

  for pt = 1:(length(x)-1)
    pslope = yslope;                  % X-axis score
    pintcp = y(pt) - pslope * x(pt);
    xp = (pintcp - xintcp)/(xslope - pslope);
    yp = xslope * xp + xintcp;
    xvalue = xscale * eucl([xlowx xlowy;xp yp]);
    if (xp < xlowx)
      xvalue = -xvalue;
    end;
    xvalue = xvalue + xlow;
    if (xlog)
      xvalue = exp(xvalue);
    end;

    pslope = xslope;                  % Y-axis score
    pintcp = y(pt) - pslope * x(pt);
    xp = (pintcp - yintcp)/(yslope - pslope);
    yp = yslope * xp + yintcp;
    yvalue = yscale * eucl([ylowx ylowy;xp yp]);
    if (yp < ylowx)
      yvalue = -yvalue;
    end;
    yvalue = yvalue + ylow;
    if (xlog)
      xvalue = exp(xvalue);
    end;

    x(pt) = xvalue;
    y(pt) = yvalue;
  end;

  x = x(1:length(x)-1);                        % Reduce to points of interest
  y = -y(1:length(y)-1);

  return;
