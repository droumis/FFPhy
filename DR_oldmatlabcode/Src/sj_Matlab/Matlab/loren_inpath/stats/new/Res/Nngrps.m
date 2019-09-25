% NNGRPS: Finds euclidean nearest-neighbor (single-linkage) distances among groups.
%
%     [dist,distmat] = nngrps(g,X)
%
%         g = vector (length n) of group identifiers for k groups.
%         X = [n x p] matrix of coordinates.
%         --------------------------------------------------------
%         dist =    [k x k] distance matrix for k groups.
%         distmat = [k(k-1)/2 x 5] matrix of pairwise distances among groups:
%                     col 1: identifier of first group;
%                         2: identifier of second group;
%                         3: index of individual within first group;
%                         4: index of individual within second group;
%                         5: distance between first and second groups.
%

% RE Strauss, 1/2/97

function [dist,distmat] = nngrps(g,X)
  grpid = uniquef(g);                 % Unique group identifiers
  ngrps = length(grpid);              % Number of groups

  dist =    zeros(ngrps,ngrps);       % Allocate output matrices
  distmat = zeros(ngrps*(ngrps-1)/2,5);
  
  k = 0;
  for i = 1:(ngrps-1)                 % Cycle thru pairs of groups
    indx1 = find(g==grpid(i));          % Isolate first group
    x1 = X(indx1,:);

    for j = (i+1):ngrps
      indx2 = find(g==grpid(j));        % Isolate second group
      x2 = X(indx2,:);

      d = eucl(x1,x2);                  % Euclidean distances between groups
      dmin = min(min(d));               % Min (nearest-neighbor) distance
      [imin,jmin] = find(d==dmin);      % Individuals having the min distance
      imin = imin(1);
      jmin = jmin(1);

      dist(i,j) = dmin;
      dist(j,i) = dmin;

      k = k+1;
      distmat(k,:) = [grpid(i) grpid(j) indx1(imin) indx2(jmin) dmin];
    end;
  end;

  return;