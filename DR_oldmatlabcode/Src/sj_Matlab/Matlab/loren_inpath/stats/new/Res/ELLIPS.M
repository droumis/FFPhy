% ELLIPS:  Computes a set of rectangular coordinates describing an ellipse
%           centered on (h,k) with major-axis length 2a and minor-axis length 2b,
%           rotated by an angle theta.
%
%       Usage: [x,y] = ellips(a,b,{h},{k},{theta},{n})
%
%           a -     half-length of major axis.
%           b -     half-length of minor axis.
%           h,k -   coordinates of center of ellipse [default = 0].
%           theta - angle of counterclockwise rotation, in radians [default = 0].
%           n -     even number of point coordinates to be returned 
%                     [default = 146].
%

% RE Strauss, 9/1/95
%   9/20/99 - update handling of null input arguments.

function [x,y] = ellips(a,b,h,k,theta,n)
  if (nargin < 3) h = []; end;
  if (nargin < 4) k = []; end;
  if (nargin < 5) theta = []; end;
  if (nargin < 6) n = []; end;

  if (isempty(h))
    h = 0;
  end;
  if (isempty(k))
    k = 0;
  end;
  if (isempty(theta))
    theta = 0;
  end;
  if (isempty(n))
    n = 146;
    N = 72;
  else
    N = floor(n/2)-1;
  end;

  x1 = -a + (a/(N/2))*(0:N);
  x2 = x1((N+1):-1:1);
  x1s = x1 .* x1;
  x2s = x2 .* x2;
  x = [x1 x2]';

  as = a * a;
  y1 = +b * abs(sqrt(1-x1s/as));
  y2 = -b * abs(sqrt(1-x2s/as));
  y = [y1 y2]';

  if (theta ~= 0)
    r = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    coord = [x y]*r;
    x = coord(:,1);
    y = coord(:,2);
  end;

  x = x+h;
  y = y+k;

  return;
