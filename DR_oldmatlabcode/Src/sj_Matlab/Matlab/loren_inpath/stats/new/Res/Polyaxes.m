% POLYAXES: Finds the anisotropy (sqrt ratio of eigenvalues) and direction of a convex 
%           polygon, based on Harvey (1981, Computers & Geosciences 7:387-392).
%
%     Syntax: [anisotropy,azimuth,eigval,axislen,cross] = polyaxes(poly,doplot)
%
%           poly =   [npt x 2] matrix of point coordinates of polygon vertices.
%           doplot = boolean flag specifying a plot of polygon and 
%                      axes [default = 0].
%           -------------------------------------------------------------------
%           anisotropy = sqrt(e1/e2), where e1 & e2 are the normalized 1st & 2nd  
%                          eigenvalues of a corresponding eigenellipse.
%           azimuth =    direction of the 1st eigenvector, counterclockwise 
%                          from horizontal, in degrees (range 0-180).  
%                          Returns NaN for singular polygon (shape=1).
%           eigval =     [2 x 1] vector of normalized eigenvalues, summing 
%                          to unity; thus they are proportions of variance 
%                          accounted for by the 1st & 2nd eigenvectors.
%           axislen =   [2 x 1] vector of lengths of major and minor axes.
%                          Product of axis lengths is scaled to total area.  
%           cross =      [4 x 2] matrix of coordinates of endpoints of cross 
%                          formed by major and minor axes, centered on 
%                          centroid.  First two rows are major axis, 
%                          second two rows are minor axis.

% RE Strauss, 4/29/97
%   9/3/99 -  changed plot colors for Matlab v5.
%   9/15/99 - added 'crosslen' to output arguments; changed name from polystrn().

function [shape,azimuth,eigval,axislen,cross] = polyaxes(poly,doplot)
  if (nargin < 2) doplot = []; end;

  if (isempty(doplot))
    doplot = 0;
  end;

  get_cross = 0;
  if (nargout > 3 | doplot)
    get_cross = 1;
  end;

  [npt,c] = size(poly);
  if (c~=2)
    if (npt==2)
      poly = poly';
      [npt,c] = size(poly);
    else
      error('  Input coordinate matrix is wrong size');
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
  
  M = (1/(2*npt))*[x2 xy; xy y2];     % Inertia matrix of 2nd moments

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

  % Find endpoints of major and minor axes for plotting.  Cross is 
  % centered on centroid.  Product of axes is scaled to total area of polygon.

  cross = [];
  if (get_cross & finite(azimuth))
    scale = 0.15;
    cross = zeros(4,2);

    A = scale * totarea;              % Distances of endpoints from centroid
    r = shape * shape;
    dist1 = sqrt(A*r) / 2;
    dist2 = A / (2*dist1) / 2;

    cross(1,1) = origcentr(1) + dist1 * cos(d1);     % Endpoints of major axis
    cross(1,2) = origcentr(2) + dist1 * sin(d1);
    cross(2,1) = origcentr(1) + dist1 * cos(d1+pi);
    cross(2,2) = origcentr(2) + dist1 * sin(d1+pi);

    cross(3,1) = origcentr(1) + dist2 * cos(d2);     % Endpoints of minor axis
    cross(3,2) = origcentr(2) + dist2 * sin(d2);
    cross(4,1) = origcentr(1) + dist2 * cos(d2+pi);
    cross(4,2) = origcentr(2) + dist2 * sin(d2+pi);

    len1 = sqrt(totarea*r) / 2;       % Lengths of axes
    len2 = totarea / (2*len1) / 2;
    axislen = [len1; len2];  
  end;

  if (doplot)
    plot(origpoly(:,1),origpoly(:,2),'ok',origpoly(:,1),origpoly(:,2),'k');
    hold on;
    if (~isempty(cross))
      plot(cross(1:2,1),cross(1:2,2),'k',cross(3:4,1),cross(3:4,2),'k');
    end;
    putbnd(origpoly(:,1),origpoly(:,2));
    axis('square');
    axis('equal');
    set(gca,'Xtick',[]);              % Suppress x-axis labels and tick marks
    set(gca,'Xcolor','k');            % Make x-axes invisible
    set(gca,'Ytick',[]);              % Suppress y-axis labels and tick marks
    set(gca,'Ycolor','k');            % Make y-axes invisible
    hold off;
    
  end;

  return;

