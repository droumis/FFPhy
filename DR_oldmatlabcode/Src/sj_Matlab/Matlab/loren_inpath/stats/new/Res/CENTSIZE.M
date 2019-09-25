% CENTSIZE: Centroid size of a point configuration.
%
%       Usage: csize = centsize(coord)
%
%           coord - [n x p] matrix of point coordinates
%           csize - centroid size: root-mean-square distance of points
%                     from their centroid
%

% RE Strauss, 8/2/95

function csize = centsize(coord)
  [n,p] = size(coord);
  csize = sqrt(mean(sum(((coord-ones(n,1)*mean(coord)).^2)')));

  return;
