% NNDIST: Finds nearest-neighbor distance for each point.
%
%     Usage: [nnd,nn] = nndist(crds,{getall},{tol})
%
%         crds =    [n x p] matrix of point coordinates in P dimensions.
%         getall =  optional boolean flag indicating, if true, that all points 
%                     having identical nn distances to a given point are to be 
%                     returned; if false, only the first in the sequence 1:n 
%                     is returned.
%         tol =     optional tolerance for judging identical distances 
%                     [default = 1e-6 x min range of points].
%         -----------------------------------------------------------------------
%         nnd =     [n x 1] vector of nearest-neighbor distances.
%         nn =      corresponding [n x 1] vector of indices of nearest neighbors, 
%                     or [n x m] matrix of indices where m is the maximum number 
%                     identical nn distances encountered for any point 
%                     of (=1 if all distances are unique).
%

% RE Strauss, 6/6/00
%   9/27/00 - allow for identical nn distances per point.

function [nnd,nn] = nndist(crds,getall,tol)
  if (nargin < 2) getall = []; end;
  if (nargin < 3) tol = []; end;

  if (isempty(getall))
    getall = 0;
  end;
  if (isempty(tol) & getall)
    tol = 1e-6 * min(range(crds));
  end;

  [n,p] = size(crds);

  dist = eucl(crds);                    % Matrix of pairwise distances
  maxdist = max(max(dist));
  for i = 1:n                           % Replace zeros on diagonal with highvals
    dist(i,i) = 2*maxdist;
  end;

  [nnd,nn] = min(dist);                 % Find min interpoint distances
  nnd = nnd';                           % Transform to col vectors
  nn = nn';   

  if (getall)
    for i = 1:n                           % Check for non-unique dists
      d = abs(dist(i,:) - nnd(i));
      j = find(d <= tol);
      lenj = length(j);
      if (lenj > 1)
        di = dist(i,j);
        nnd(i) = mean(di);
        [nn,j] = padcols(nn,j);
        nn(i,:) = j;
      end;      
    end;
  end;

  return;
