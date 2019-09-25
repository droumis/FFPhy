% MakePolygon: Given a set of 2D points, arranges them into a counterclockwise path
%              (with respect to the centroid) describing a closed polygon.
%
%     Usage: poly = makepolygon(pts,doplot)
%
%         pts =     [n x 2] matrix of point coordinates.
%         --------------------------------------------------------------------
%         poly =    [n+1 x 2] matrix of vertices of a closed polygon.
%         doplot =  optional boolean variable indicating, if true, that a plot
%                     of the polygon is to be produced [default = 0].
%

% RE Strauss, 5/21/03

function poly = makepolygon(pts,doplot)
  if (~nargin) help makepolygon; return; end;
  
  if (nargin < 2) doplot = []; end;
  if (isempty(doplot)) doplot = 0; end;
  
  centroid = mean(pts);
  theta = anglerotation(centroid,pts(1,:),pts);   % Angles from first point in list
  [theta,poly] = sortmat(theta,pts);              % Sort points by theta
  poly = [poly; poly(1,:)];                       % Repeat first point at end
  
  if (doplot)
    figure;
    plot(poly(:,1),poly(:,2),'k',poly(:,1),poly(:,2),'ko',...
         poly(1,1),poly(1,2),'k*',centroid(1),centroid(2),'k+');
    putbnds(poly);
  end;
  
  return;
  