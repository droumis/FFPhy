% RANDCNVX: Generates uniform-random points within a given CONVEX polygon.  
%           Assumes convexity, doesn't check for it.
%
%     Syntax: pt = randcnvx(poly,npts,doplot)
%
%           poly -   [m x 2] matrix of coordinates specifying a polygon.
%           npts -   number of random points to be generated [default=1].
%           doplot - boolean value indicating whether (1,=TRUE) or 
%                     not (0,=FALSE) to produce a plot of the polygon and 
%                     points [default=0].
%           -----------------------------------------------------------
%           pt -     [npts x 2] matrix of random point coordinates.
%

function pt = randcnvx(poly,npts,doplot)
  if (nargin < 2)                         % Default number of points
    npts = 1;
  end;
  if (nargin < 3)
    doplot = 0;
  end;

  [m,p] = size(poly);
  if (p~=2)                               % Check dimensions
    error('  Procedure implemented only for 2 dimension');
  end;

  R = rand(npts,2);                       % Uniform-random points

  if (poly(m,:)~=poly(1,:))               % Close polygon if necessary
    poly = [poly; poly(1,:)];
    m = m+1;
  end;

  p = 1:(m-1);
  q = 2:m;
  edge = [poly(p,2) - poly(q,2), ...      % Polygon edges
          poly(q,1) - poly(p,1), ...
          poly(p,1).*poly(q,2) - poly(p,2).*poly(q,1)];

  poly = poly(p,:);                       % Open polygon
  m = m-1;

  centr = mean(poly);                     % Centroid

  crds = (poly - (ones(m,1)*centr))';     % Center points on centroid
  adev = angledev(poly(1,:),centr,poly,1);  % Cum spoke angles about centroid
  cumangle = [adev(2:m); 2*pi];
  angle = [cumangle(1);cumangle(2:m)-cumangle(1:(m-1))]; % Diffs

  randangle = (2*pi*R(:,1));              % Random angles about circle
  ra =  randangle* ones(1,m);             % Propagate columns
  ang = ones(npts,1) * adev(1:m)';        % Propogate columns
  dev = ra - ang;                         % Angular deviations from each spoke
  dev = dev(:);                           % Reassign neg angles to high value
  indx = find(dev<0);
  dev(indx) = 10*ones(length(indx),1);
  dev = reshape(dev,npts,m)';
  [dev,t] = min(dev);                     % Random angular deviations past 't' spokes
  dev = dev';
  t = t';

  theta = angledev(centr+[1 0],centr,poly(1,:),1);  % Ang deviation from abscissa
  adjrandangle = randangle + theta*ones(npts,1);  % Adjusted random angles

  c = ones(npts,1)*centr;                 % Propagated centroid
  endpt = c + [cos(adjrandangle), sin(adjrandangle)];  % Unit radius endpt
  radius = [c(:,2) - endpt(:,2), ...      % Randomized radii
            endpt(:,1) - c(:,1), ...
            c(:,1).*endpt(:,2) - c(:,2).*endpt(:,1)];
  bndry = edge(t,:);                      % Edges for randomized points

  dlm = radius(:,1).*bndry(:,2) - radius(:,2).*bndry(:,1);
  intsect = [(radius(:,2).*bndry(:,3) - radius(:,3).*bndry(:,2))./dlm, ...
             (radius(:,3).*bndry(:,1) - radius(:,1).*bndry(:,3))./dlm];  % Intsect w/ bndry
  radius = sqrt((c(:,1)-intsect(:,1)).^2 + (c(:,2)-intsect(:,2)).^2); % Length radius
  randlen = radius .* sqrt(R(:,2));       % Random length along radius

  pt = c + [randlen.*cos(adjrandangle), randlen.*sin(adjrandangle)];

  if (doplot)
    clf;
    hold on;
    plot(poly(:,1),poly(:,2),'y-',pt(:,1),pt(:,2),'ro');
    plot(poly(:,1),poly(:,2),'yo');
    plot([poly(m,1),poly(1,1)],[poly(m,2),poly(1,2)],'y-');
    axis('equal');
    hold off;
  end;

  return;

