% TRIANGPT: Triangulates the coordinates of a point (p3) in relation to two others (p1, p2),
%           given the three interpoint distances.
%           
%
%     Syntax: P3 = triangpt(crds,dists)
%
%           crds =  [2 x 2] matrix of coordinates of points p1 and p2:
%                       x1,y1
%                       x2,y2
%           dists = vector of the interpoint distances from the existing 
%                     points p1 and p2 to the new point p3:
%                       1: d13 =  distance from p1 to p3.
%                       2: d23 =  distance from p2 to p3.
%           -----------------------------------------------------------------------
%           P3 =    [2 x 2] matrix of the new point p3 and its reflection about the 
%                     line p1-p2.  The first point is to the left of the ray p1-p2,
%                     and its reflection is to the right.
%                     
%

% RE Strauss, 3/5/96
%   5/2/03 -  correct the final rotation procedure; change output format.
%   5/11/03 - delete warning message.

function P3 = triangpt(crds,dists)
  d12 = eucl(crds);
  d13 = dists(1);
  d23 = dists(2);

  % p1 is placed at (0,0), p2 at (d12,0)

  [newpts,theta] = regrot(1,2,crds);

  x1 = newpts(1,1);
  y1 = newpts(1,2);
  x2 = newpts(2,1);
  y2 = newpts(2,2);

  % p3 lies at the intersection of a circle of radius d13 about p1 --
  %       x^2 + y^2 = (d13)^2
  % -- and the circle of radius d23 about point 2 --
  %       (x-x2)^2 + y^2 = (d23)^2.  
  % The difference between these equations is the equation of a line --
  %       2*x2 * x = (d13)^2 - (d23)^2 + x2^2
  % -- the so-called radical axis of the two circles, which includes
  % their intersections.  Solve for x as the x-coordinate of the new 
  % point 3.

  x3 = (d13^2 - d23^2 + x2^2) / (2*x2);

  % If the x-coordinate is farther from the origin than the radii
  % of the circles about p1 or p2, the construction is impossible.
  
  if (abs(x3) >= d13 | abs(x3-x2) >= d23)
    x = [NaN;NaN];
    y = [NaN;NaN];
    P3 = NaN;
%     disp('  TRIANGPT Warning: triangulation failed');
    return;
  end;

  % Find the y-coordinate of point 3 from the above equation --
  %       x^2 + y^2 = (d13)^2

  y3 = sqrt(d13^2 - x3^2);

  % Rotate and translate p3 and its reflection into the original 
  % coordinate system.
  
  x3 = [x3;x3];
  y3 = [y3;-y3];
  P3 = rotate([x3 y3],-theta,[0 0]) + crds([1,1],:);

%   figure;
%   plot(crds(:,1),crds(:,2),'b',P3(:,1),P3(:,2),'bo');
%   putbnds([crds(:,1);P3(:,1)],[crds(:,2);P3(:,2)]);
%   axis equal;

  return;
