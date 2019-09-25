% NOT COMPLETE

% CURVATURE:  Calculates the curvature (reciprocal of radius) at each point for 
%             a closely spaced series of points describing a 2D boundary.  
%             Estimated at each vertex from the circle that is tangent to the 
%             midpoints of the two edges.  If the coordinates do not describe a 
%             closed polygon, the curvature of the open endpoints is defined to 
%             be zero.
%
%     Usage: curv = curvature(crds)
%
%             crds = [n x 2] matrix of point coordinates.
%             -----------------------------------------------------------
%             curv = [n x 1] vector of corresponding curvature estimates.
%

% RE Strauss, 3/23/00

function curv = curvature(crds)
  [n,p] = size(crds);
  if (p ~= 2)
    error('  CURVATURE: coordinates must be two-dimensional.');
  end;

  curv = zeros(n,1);

  for i = 2:(n-1)                         % Curvature at each vertex
i
    [p1,p2] = extrcols(crds(i-1,:));
    [q1,q2] = extrcols(crds(i,:));
    [r1,r2] = extrcols(crds(i+1,:));
points = [p1 p2; q1 q2; r1 r2]

    line1 = [p2-q2, q1-p1, p1*q2-q1*p2];
    line2 = [q2-r2, r1-q1, q1*r2-r1*q2];
    lenline1 = eucl([p1 p2],[q1 q2]);
    lenline2 = eucl([q1 q2],[r1 r2]);
%lines = [line1; line2]

    mid1 = [(p1+q1)/2, (p2+q2)/2];
    mid2 = [(q1+r1)/2, (q2+r2)/2];
midpoints = [mid1; mid2]

    perp1 = [-line1(2), line1(1), line1(2)*mid1(1)-line1(1)*mid1(2)];
    perp2 = [-line2(2), line2(1), line2(2)*mid2(1)-line2(1)*mid2(2)];
%perps = [perp1; perp2]

    dlm = perp1(1)*perp2(2) - perp1(2)*perp2(1);
dlm
    if (dlm > 1e-6)
      intsct = [(perp1(2)*perp2(3)-perp1(3)*perp2(2))/dlm, ...
                (perp1(3)*perp2(1)-perp1(1)*perp2(3))/dlm];
intsct
      radius = eucl(mid1,intsct); %sqrt((mid1(1)-intsct(1)).^2+(mid1(2)-intsct(2)).^2);
rad1 = radius;
rad2 = sqrt((mid2(1)-intsct(1)).^2+(mid2(2)-intsct(2)).^2);
rad2 = eucl(mid2,intsct);
rad = [rad1 rad2]
      curv(i) = 1./radius;
    end;
  end;

  if (eucl(crds(1,:),crds(n,:)) < 1e-6)   % Closed polygon
    [p1,p2] = extrcols(crds(n-1,:));
    [q1,q2] = extrcols(crds(1,:));
    [r1,r2] = extrcols(crds(2,:));

    line1 = [p2-q2, q1-p1, p1*q2-q1*p2];
    line2 = [q2-r2, r1-q1, q1*r2-r1*q2];

    mid1 = [(p1+q1)/2, (p2+q2)/2];
    mid2 = [(q1+r1)/2, (q2+r2)/2];

    perp1 = [-line1(2), line1(1), line1(2)*mid1(1)-line1(1)*mid1(2)];
    perp2 = [-line2(2), line2(1), line2(2)*mid2(1)-line2(1)*mid2(2)];

    dlm = line1(1)*line2(2) - line1(2)*line2(1);
    if (dlm > 1e-6)
      intsct = [(line1(2)*line2(3)-line1(3)*line2(2))/dlm, ...
                (line1(3)*line2(1)-line1(1)*line2(3))/dlm];
      radius = sqrt((mid1(1)-intsct(1)).^2+(mid1(2)-intsct(2)).^2);
      curv(1) = 1./radius;
      curv(n) = curv(1);
    end;
  end;

  return;
