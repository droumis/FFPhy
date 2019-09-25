% CIRCFITF: Objective function for circfit, which fits a circle to a set of 
%           2-dimensional points by finding the center point having minimum variance 
%           in Euclidean distance to all points.
%
%     Usage: vardist = circfitf(center,crds)
%
%           center = [1 x 2] row vector of point coordinates of fitted center.
%           crds =   [n x 2] matrix of point coordinates.
%           ------------------------------------------------------------------
%           var =    variance of Euclidean distances.
%

function vardist = circfitf(center,crds)
  vardist = var(eucl(crds,center));
  return;


