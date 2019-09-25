% POLYAREA: Area, perimeter, and centroid of a 2-dimensional polygon.
%
%     Syntax: [area,perim,centroid] = polyarea(crds)
%
%           crds =     [n x 2] matrix of coordinates of the polygon vertices.
%           -----------------------------------------------------------------
%           area =     polygon area.
%           perim =    polygon perimeter.
%           centroid = [1 x 2] vector of centroid coordinates.
%

% Harvey, P.K. 1981. A simple algorithm for the unique characterization of 
%   convex polygons.  Comput. Geosci. 7:387-392.

% RE Strauss, 5/7/97
%   2/3/00 -  corrected error for non-convex polygons.
%   3/19/00 - corrected error for centroid when vertices are clockwise.
%   3/14/02 - return NaN's for if too few coordinates are passed.
%   5/6/03 -  return reasonable values for two points.

function [area,perim,centroid] = polyarea(crds)
  [n,p] = size(crds);
  if (p~=2)                     % Check matrix size
    if (n==2)
      crds = crds';
      [n,p] = size(crds);
    else
      error('  POLYAREA: input matrix is wrong size');
    end;
  end;
  
  if (n<2)
    area = NaN;
    perim = NaN;
    centroid = NaN;
    return;
  elseif (n==2)
    area = 0;
    perim = 2*eucl(crds);
    centroid = mean(crds);
  end;

  if (crds(1,:)==crds(n,:))     % If polygon is closed, open it
    n = n-1;
  end;

  t = 0;
  tx = 0;
  ty = 0;

  x = crds(1,1);
  y = crds(1,2);

  x2 = crds(2,1);
  y2 = crds(2,2);

  perim = eucl([x y],[x2 y2]) + eucl([x y],crds(n,:));

  for i = 2:(n-1)
    x1 = x2;
    x2 = crds(i+1,1);
    y1 = y2;
    y2 = crds(i+1,2);

    perim = perim + eucl([x1 y1],[x2 y2]);

    xmid = (x1+x2)/2;
    ymid = (y1+y2)/2;
    xd = (xmid-x)*2/3 + x;
    yd = (ymid-y)*2/3 + y;
    area = (x*y1 - x*y2 + x1*y2 - x1*y + x2*y - x2*y1)/2;

    t = t+area;
    tx = tx + area*xd;
    ty = ty + area*yd;
  end;

  if (abs(t)>eps)
    centroid = [tx/t ty/t];
  else
    centroid = mean(crds);
  end;
  area = abs(t);

  return;
