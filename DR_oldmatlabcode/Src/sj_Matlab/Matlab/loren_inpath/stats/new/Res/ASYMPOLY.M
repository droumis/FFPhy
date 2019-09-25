% ASYMPOLY: Given a polygon described by a configuration of points, measures
%           the minimum bilateral asymmetry from an initial partitioning by
%           the major axis, optimized so that the reflection of the polygon about the 
%           axis is optimally superpositioned on the original.  
%             Two asymmetry measures are provided: difference in half-area and 
%           difference in half-perimeter, where each difference is the 'first-half'
%           minus the 'second-half' (in terms of the sequence of points).
%
%     Usage: [asym_area,asym_perim] = asympoly(pt,{doplot})
%
%           pt =         [p x 2] matrix of point coordinates.
%           doplot =     optional boolean flag indicating that plots are to be produced 
%                          [default = 0].
%           -------------------------------------------
%           asym_area =  [2 x 1] vector of difference in half-areas: absolute, and 
%                          relative to total area.
%           asym_perim = correslponding [2 x 1] vector of differences in half-perimeters.
%
%

% RE Strauss, 6/26/98
%   8/19/99 - change plot colors for Matlab v5.
%   1/4/00 -  changed fminu() to fmins().

function [asym_area,asym_perim] = asympoly(pts,doplot)
  if (nargin < 2) doplot = []; end;

  if (isempty(doplot))
    doplot = 0;
  end;

  [npts,c] = size(pts);

  b = majaxis(pts(:,1),pts(:,2));     % Major axis is initial guess at midsaggital line
  init_line = b(1,:);

  if (abs(init_line(1)) > 5);         % If slope too steep
    pts = rotate(pts,-pi/2);          %   rotate by 90 degrees
    b = majaxis(pts(:,1),pts(:,2));   % New major axis
    init_line = b(1,:);
  end;

  line = fmins('refldiff',init_line,[],[],pts);   % Optimize position of line

  if (doplot)
    b1 = line(1);
    b0 = line(2);
    xmin = min(pts(:,1));
    xmax = max(pts(:,1));
    int1 = [xmin b0+b1*xmin];
    int2 = [xmax b0+b1*xmax];
    l = [int1; int2];

    close all;
    plot(pts(:,1),pts(:,2),'yo',pts(:,1),pts(:,2),'b-',l(:,1),l(:,2),'k-');
    putbnd(pts(:,1),pts(:,2));
    axis('equal');
    reflplot(line,pts);
  end;

  asym_area = zeros(2,1);
  asym_perim = zeros(2,1);

  [asym_area(1),asym_perim(1)] = reflstat(line,pts);    % Calculate lateral asymmetries
  [area,perim] = polyarea(pts);                         % Total area and perimeter
  asym_area(2) = asym_area(1) / area;
  asym_perim(2) = asym_perim(1) / perim;

  return;

