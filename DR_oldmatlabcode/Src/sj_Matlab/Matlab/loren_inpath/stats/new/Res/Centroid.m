% CENTROID: Finds a robust 50%ile trimmed estimate of the centroid, and a 0%ile 
%           estimate of the medioid, of two-dimensional point coordinates 
%           from the central half of the data, as determined by convex-hull 
%           peeling (Riani Zani & Corbellini 1997).
%
%           For the usual parametric centroid, use polyarea().
%
%     Usage: [cent,med] = centroid(crds)
%
%           crds = [n x 2] matrix of point coordinates.
%           -------------------------------------------------------
%           cent = [1 x 2] vector of coordinates of the centroid.
%           med =  [1 x 2] vector of coordinates of the medioid.
%

% Riani, M., S. Zani & A. Corbellini. 1997. Robust bivariate boxplots and
%   visualization of multivariate data.  In: I. Balderjahn, R. Mathar &
%   M. Schader (eds.), Classification, Data Analysis, and Data Highways,
%   pp. 93-100.  Springer-Verlag.

% RE Strauss, 1/5/99
%   7/27/00 - allow for observed min pctile >50.

function [cent,med] = centroid(crds)
  get_median = 0;
  if (nargout > 1)
    get_median = 1;
  end;

  [n,p] = size(crds);
  if (p~=2)
    error('  CENTROID: coordinates must be two-dimensional');
  end;

  pctile = hullpeel(crds);              % Find point percentiles by hull-peeling
  i = find(pctile<=max([50,min(pctile)])); % Isolate central half of points
  if (length(i)>1)
    cent = mean(crds(i,:));             % Centroid
  else
    cent = crds(i,:);
  end;

  if (get_median)
    minpctile = min(pctile);            % Find central hull
    i = find(pctile==minpctile);        % Find central points
    if (length(i)>1)
      med = mean(crds(i,:));            % Medioid is mean of most central points
    else
      med = crds(i,:);
    end;
  end;

  return;
