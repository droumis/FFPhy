% RANDPATH: Generates points at uniform-random distances along a specified 
%           line-segment path.
%
%     Usage: [pt,dist] = randpath(path,{npts})
%
%           path = [m x 2] matrix of coordinates specifying a path.
%           npts = optional number of random points to be generated [default=1].
%           --------------------------------------------------------------------
%           pt =   [npts x 2] matrix of random point coordinates.
%           dist = [npts x 1] vector of corresponding distances along path.
%

% RE Strauss, 9/21/95

function [pt,dist] = randpath(path,npts)
  if (nargin < 2) npts = []; end;

  if (isempty(npts))                  % Default number of points
    npts = 1;
  end;

  [m,p] = size(path);
  pt = zeros(npts,2);

  pathlen = 0;                        % Total path length
  for mi = 1:(m-1)
    lineseg = path((mi:mi+1),:);
    pathlen = pathlen + eucl(lineseg);
  end;

  dist = rand(npts,1)*pathlen;        % Random dists along path
  for ni = 1:npts                     % Generate specified number of pts
    d = 0;
    mi = 0;
    while (d < dist(ni))                % Sum length along path segments
      mi = mi+1;                        %   till go too far
      lineseg = path((mi:mi+1),:);
      delta = eucl(lineseg);
      d = d + delta;
    end;

    d = d - delta;                      % Back up
    diff = dist(ni)-d;                  % Diff between what have and what need
    prop = diff / delta;                % Proportion of way to go to next pt
    pt(ni,:) = path(mi,:) + prop*(path(mi+1,:)-path(mi,:));
  end;

  return;

