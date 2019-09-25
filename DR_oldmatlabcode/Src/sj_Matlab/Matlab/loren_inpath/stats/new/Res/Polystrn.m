% POLYSTRN: Finds the anisotropy (sqrt ratio of eigenvalues) and direction of a 
%           convex polygon.
%
%     Syntax: [anisotropy,azimuth,eigval,cross] = polystrn(poly,doplot)
%
%           poly =   [npt x 2] matrix of point coordinates of polygon vertices.
%           doplot = boolean flag specifying a plot of polygon and 
%                      dilatations [default = 0].
%           -------------------------------------------------------------------
%           anisotropy = sqrt(e1/e2), where e1 & e2 are the normalized 1st & 2nd  
%                          eigenvalues of a corresponding eigenellipse.
%           azimuth =    direction of the 1st eigenvector, counterclockwise 
%                          from horizontal, in degrees (range 0-180).  
%                          Returns NaN for singular polygon (shape=1).
%           eigval =     [2 x 1] vector of normalized eigenvalues, summing 
%                          to unity; thus they are proportions of variance 
%                          accounted for by the 1st & 2nd eigenvectors.
%           cross =      [4 x 2] matrix of coordinates of endpoints of cross 
%                          formed by major and minor dilatations, centered on 
%                          centroid.  First two rows are major dilatation, second 
%                          two rows are minor dilatation.
%

% Harvey, 1981, Computers & Geosciences 7:387-392.

% RE Strauss, 4/29/97
%   5/21/03 - modified documentation, initial 'nargin' statements, error messages;
%             modified figure from Matlab4 to Matlab6 conventions.

function [shape,azimuth,eigval,cross] = polystrn(poly,doplot)
  if (~nargin) help polystrn; return; end;
  
  if (nargin < 2) doplot = []; end;
  
  if (isempty(doplot)) doplot = 0; end;

  if (nargout > 3 | doplot)
    get_cross = 1;
  else
    get_cross = 0;
  end;

  [npt,c] = size(poly);
  if (c~=2)
    if (npt==2)
      poly = poly';
      [npt,c] = size(poly);
    else
      error('  POLYSTRN: Input coordinate matrix is wrong size');
    end;
  end;

  if (poly(npt,:)~=poly(1,:))       % If first vertex not repeated at end,
    poly = [poly; poly(1,:)];       %   add onto end
  else                              % else
    npt = npt-1;                    %   decrement vertex count
  end;

  [totarea,perim,centr] = polyarea(poly);  % Centroid of convex polygon
  origpoly = poly;
  origcentr = centr;

  poly = poly - ones(npt+1,1)*centr;  % Zero-center
  centr = [0 0];

  area = zeros(npt,1);              % Allocate matrices
  fP = zeros(2*npt,1);
  P = zeros(2*npt,2);
  
  for t = 1:npt                     % Cycle thru triangles
    area(t) = polyarea([poly(t:(t+1),:);centr]);    % Area

    P(t*2-1,:) = mean([poly(t,:);centr]); % Interior
    if (t>1)
      fP(t*2-1) = sum(area((t-1):t));
    end;

    P(t*2,:) = mean(poly(t:(t+1),:));       % Boundary
    fP(t*2) = area(t);
  end;
  fP(1) = area(1)+area(npt);
  fP = fP/3;

  x2 =  sum(fP.*P(:,1).*P(:,1));
  y2 =  sum(fP.*P(:,2).*P(:,2));
  xy = -sum(fP.*P(:,1).*P(:,2));
  
  M = (1/(2*npt))*[x2 xy; xy y2];  % Inertia matrix of 2nd moments

  [evects,evals] = eigen(M);
  eigval = evals/sum(evals);
  shape = sqrt(evals(1)/evals(2));    % Shape statistic

  if (abs(shape-1) > 1e-4)            % If not singular,
    d1 = atan(evects(1,2)/evects(2,2)); % Counterclockwise azimuths from horizontal
    d2 = atan(evects(1,1)/evects(2,1));
    if (d1 < 0)
      d1 = d1 + pi;
    end;
    if (d2 < 0)
      d2 = d2 + pi;
    end;

    deg_per_rad = 360/(2*pi);           
    azimuth = d1 * deg_per_rad;
  else
    azimuth = NaN;
  end;

  % Find endpoints of major and minor dilatations for plotting.  Cross is 
  % centered on centroid.  Product of dilitations is scaled to total area of polygon.

  cross = [];
  if (get_cross & finite(azimuth))
    scale = 0.15;
    cross = zeros(4,2);

    A = scale * totarea;              % Distances of endpoints from centroid
    r = shape * shape;
    dist1 = sqrt(A*r) / 2;
    dist2 = A / (2*dist1) / 2;

    cross(1,1) = origcentr(1) + dist1 * cos(d1);     % Endpoints of major dilatation
    cross(1,2) = origcentr(2) + dist1 * sin(d1);
    cross(2,1) = origcentr(1) + dist1 * cos(d1+pi);
    cross(2,2) = origcentr(2) + dist1 * sin(d1+pi);

    cross(3,1) = origcentr(1) + dist2 * cos(d2);     % Endpoints of minor dilatation
    cross(3,2) = origcentr(2) + dist2 * sin(d2);
    cross(4,1) = origcentr(1) + dist2 * cos(d2+pi);
    cross(4,2) = origcentr(2) + dist2 * sin(d2+pi);
  end;

  if (doplot)
    plot(origpoly(:,1),origpoly(:,2),'ok',origpoly(:,1),origpoly(:,2),'k');
    hold on;
    if (~isempty(cross))
      plot(cross(1:2,1),cross(1:2,2),'k',cross(3:4,1),cross(3:4,2),'k');
    end;
    putbnd(origpoly(:,1),origpoly(:,2));
    axis square;
    axis equal;
    axis off;
    hold off;
  end;

  return;

