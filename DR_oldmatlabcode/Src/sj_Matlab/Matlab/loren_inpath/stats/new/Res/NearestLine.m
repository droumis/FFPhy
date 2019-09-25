% NearestLine: finds the nearest line segment for one or more points, in 2D space.
%
%     Usage: nearest = nearestline(pt,segments)
%
%         pt = [n x 2] matrix of point coordinates (in p-dimensional space) for
%                 n points.
%         segments = [m+1 x 2] matrix of point coordinates specifying a path
%                 containing m line segments.
%         ------------------------------------------------------------------------
%         nearest = [n x 1] vector of integers (range 1-m) identifying the nearest
%                 line segment for each of n points.
%

% RE Strauss, 10/16/02

function nearest = nearestline(pt,segments)
  [n,p] = size(pt);
  [m1,q] = size(segments);
  
  if (p~=q)
    error('  NearestLine: points and line segments must be of equal dimension.');
  end;
  
  m = m1-1;
  
  dist = zeros(m,n);              % Distances from points to line segments
  for iseg = 1:m
    dist(iseg,:) = distptsline(segments(iseg:(iseg+1),:),pt)';
  end;
  
  [d,nearest] = min(dist);
  nearest = nearest';

  return;
  