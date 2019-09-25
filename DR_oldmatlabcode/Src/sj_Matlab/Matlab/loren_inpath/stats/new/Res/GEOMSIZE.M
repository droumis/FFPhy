% GEOMSIZE: Geometric size of a point configuration, the geometric mean 
%           of all possible interpoint distances.
%
%       Usage: gsize = geomsize(coord)
%
%           coord - [n x p] matrix of point coordinates
%           gsize - centroid size: root-mean-square distance of points
%                     from their centroid
%

% RE Strauss, 8/2/95

function gsize = geomsize(coord)
  [n,p] = size(coord);
  dist = trilow(eucl(coord));
  gsize = exp(mean(log(dist)));

  return;
