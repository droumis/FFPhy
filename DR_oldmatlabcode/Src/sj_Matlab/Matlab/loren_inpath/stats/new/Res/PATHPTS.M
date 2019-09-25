% PATHPTS:  Evenly distributes a given number of points along a path 
%           specified by point coordinates in p dimensions.
%
%     [crds,p] = pathpts(path,npts)
%
%           path = [n x 2] points defining path.
%           npts = number of points to be distributed [default = 100].
%           ----------------------------------------------------------
%           crds = coordinates of distributed points.
%           p =    cumulative distance along path for each point.
%

% RE Strauss, 8/23/98
%   3/21/00 - extend from 2 to p dimensions.

function [crds,p] = pathpts(path,npts)
  if (nargin < 2) npts = []; end;

  if (isempty(npts))
    npts = 100;
  end;

  [n,p] = size(path);

  d = [path(2:n,:) - path(1:(n-1),:)]';   % Distances between consecutive points
  d = sqrt(sum(d.^2))';

  cumpathlen = [0; cumsum(d)];            % Cum distances to points
  pathlen = sum(d);    

  crds = zeros(npts,p);
  crds(1,:) = path(1,:);
  crds(npts,:) = path(n,:);
  p = linspace(0,pathlen,npts)';

  for i = 2:npts-1
    ip = max(find(cumpathlen <= p(i)));
    dp = (p(i)-cumpathlen(ip)) ./ (cumpathlen(ip+1)-cumpathlen(ip));
    crds(i,:) = (dp.*(path(ip+1,:)-path(ip,:))) + path(ip,:);
  end;

  return;
