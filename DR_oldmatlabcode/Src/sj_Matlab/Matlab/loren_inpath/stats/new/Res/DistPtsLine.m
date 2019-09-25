% DistPtsLine: Finds distances from points to a line specified by two other points.
%
%     Usage: dist = distptsline(line,pts)
%
%         line = [2 x 2] matrix of coordinates for two points specifying the position
%                   of the line.
%         pts =  [n x 2] matrix of point coordinates.
%         ---------------------------------------------------------------------------
%         dist = [n x 1] vector of distances of points from line.
%

% RE Strauss, 3/14/02

function dist = distptsline(line,pts)
  if (size(line) ~= [2 2])
    error('  DistPtsLine: invalid coordinates for line.');
  end;
  if (size(pts,2) ~= 2)
    error('  DistPtsLine: invalid coordinates for points.');
  end;

  crds = [line; pts];
  crds = regrot(1,2,crds);
  dist = abs(crds(3:end,2));
  return;
  
  

