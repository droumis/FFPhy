% POLYGON: Finds the area, perimeter, and length of side of a regular n-sided 
%          polygon (n-gon) of given radius (center to vertex).
%
%     Usage: [area,perim,sidelen] = polygon(n,{radius},{relative})
%
%          n =        order of polygon.
%          radius =   optional radius [default=1].
%          relative = optional boolean flag indicating, if true (=1), that 
%                       area, perimeter, and side length are to be reported as 
%                       proportions of values for a circle [default=0].
%          -------------------------------------------------------------------
%          area =     polygon area.
%          perim =    polygon perimeter.
%          sidelen =  length of side of polygon.
%

function [area,perim,sidelen] = polygon(n,radius,relative)
  if (nargin < 2)
    radius = [];
  end;
  if (nargin < 3)
    relative = [];
  end;

  if (isempty(radius))
    radius = 1;
  end;
  if (isempty(relative))
    relative = 0;
  end;

  area = 0.5*n*radius*radius*sin(2*pi/n);
  perim = 2*n*radius*sin(pi/n);
  sidelen = perim/n;

  if (relative)
    A = pi*radius*radius;
    P = 2*pi*radius;
    area = area/A;
    perim = perim/P;
    sidelen = sidelen/P;
  end;

  return;
