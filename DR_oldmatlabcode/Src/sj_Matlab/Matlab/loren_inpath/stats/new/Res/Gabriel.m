% GABRIEL: Finds the Gabriel connectivity graph among points in P dimensions.
%
%     Usage: [connect,dist] = gabriel(crds,{noplot},{tol})
%
%           crds =    [n x p] matrix of point coordinates.
%           noplot =  optional boolean flag indicating that plot should be 
%                       suppressed [default = 0].  Plot of the Gabriel graph
%                       is produced only for p=2.
%           tol =     optional tolerance for squared distance from a third point 
%                       to the minimum A-B circle (i.e., that circle on whose 
%                       circumference A & B are at opposite points) [default = 1e-6].
%           -------------------------------------------------------------------------
%           connect = [n x n] boolean adjacency matrix.
%           dist =    corresponding edge lengths (Euclidean distances);
%                       non-connected edge distances are given as zero.
%

% Gabriel, KR & RR Sokal. 1969. A new statistical approach to geographic variation
%   analysis.  Syst. Zool. 18:259-278.
% Matula, DW & RR Sokal. 1980. Properties of Gabriel graphs relevant to geographic
%   variation research and the clustering of points in the plane.  Geogr. Analysis
%   12:205-222.

% RE Strauss, 1/8/99
%   9/7/99 -   changed plot colors for Matlab v5.
%   10/23/02 - remove restriction of p==2;
%              produce plot only for p==2.
%   6/4/03 -   added tolerance to the minimum A-B circle.

function [connect,dist] = gabriel(crds,noplot,tol)
  if (nargin < 2) noplot = []; end;
  if (nargin < 3) tol = []; end;

  if (isempty(noplot)) noplot = 0; end;
  if (isempty(tol))    tol = 1e-6; end;

  [n,p] = size(crds);

  if (n<2)
    error('  GABRIEL: requires at least two points ');
  end;

  connect = zeros(n,n);

  dist = eucl(crds).^2;                 % Find all squared pairwise interpoint distances
  d = dist;

  for i = 1:(n-1)                       % Cycle thru all possible pairs of points
    for j = (i+1):n
      dij = d(i,j);
      c = 1;
      for k = 1:n
        if (k~=i & k~=j)
          if (d(i,k)+d(j,k)-tol <= dij)
            c = 0;
          end;
        end;
      end;
      if (c)
        connect(i,j) = 1;
        connect(j,i) = 1;
      else
        dist(i,j) = 0;
        dist(j,i) = 0;
      end;
    end;
  end;

  if (~noplot & p==2)
    figure;
    plot(crds(:,1),crds(:,2),'ko');
    putbnds(crds(:,1),crds(:,2));
    axis('equal');
    hold on;
    for i = 1:(n-1)
      for j = 2:n
        if (connect(i,j))
          plot(crds([i j],1),crds([i j],2),'k');
        end;
      end;
    end;
    hold off;
  end;

  return;
