% QUADSPLN:  Fit a smoothed polynomial boundary to an open or closed 
%           (polygonal) 2-dimensional path.  The boundary is fitted using 
%           splined quadratic functions fitted to overlaping triplets 
%           of points (Overhauser curves).  The points returned are 
%           equally spaced along the boundary.  
%           Does not use interior center point.
%
%     Usage: scrds = quadspln(crds,{start},{doplot},{nsp})
%
%           crds =    [n x 2] set of point coordinates describing the path or 
%                       boundary.  If the first point is repeated as the last  
%                       point, the boundary is assumed to be closed.
%           start =   optional index (subscript) of starting point 
%                       [default = first point in 'crds'].
%           doplot =  optional boolean vector indicating that plot of points  
%                       and fitted boundary is to be produced [default=0].
%           nsp =    optional number of fitted points to be returned 
%                       [default=200].
%           -----------------------------------------------------------------
%           scrds =   [nsp x 2] matrix of smoothed fitted points.
%

% RE Strauss, 8/28 1998.
%   6/13/99 - add optional plot.
%   8/21/99 - changed plot colors for Matlab v5.
%   3/18/00 - changed output to single matrix rather than two vectors;
%             allowed for specification of start point.
%   3/26/00 - corrected overshot problem in neighborhoods of endpoints.

% Bowyer, A & J Woodwark. 1983. A programmer's geometry.  Butterworths, London.

function scrds = quadspln(crds,start,doplot,nsp)
  if (nargin < 2) start = []; end;
  if (nargin < 3) doplot = []; end;
  if (nargin < 4) nsp = []; end;

  [n,p] = size(crds);
  if (p ~= 2)
    error('  QUADSPLN: 2-dimensional coordinates only.');
  end;
  if (n < 3)
    error('  QUADSPLN: minimum number of points is 3.');
  end;

  if (isempty(nsp))
    nsp = 200;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;
  if (isempty(start))
    start = 1;
  end;

  if (start<1 | start>n)
    error('  QUADSPLN: start point out of range.');
  end;

  if (start > 1)                            % Shift around start point
    closed = 0;                               % Open the polygon
    if (eucl(crds([1,n],:)) < (1e-6)*max(range(crds)))
      closed = 1;
      crds(n,:) = [];
      n = n-1;
    end;

    crds = [crds(start:n,:); crds(1:start-1,:)];  % Shift

    if (closed)                               % Reclose the polygon
      crds = [crds; crds(1,:)];
      n = n+1;
    end;
  end;

  closed = 0;
  ntriplets = n-2;                          % Number of overlapping triplets of points
  if (abs(sum(crds(1,:)-crds(n,:))) < eps)  % Path closed?
    closed = 1;
    ntriplets = n-1;
    n = n+1;
    crds = [crds; crds(2,:)];
  end; 

  d = crds(1:n-1,:) - crds(2:n,:);
  d = sqrt(d(:,1).*d(:,1)+d(:,2).*d(:,2));  % Distances between consecutive points
  cumpathlen = [0; cumsum(d)];              % Cum distances to points
  pathlen = sum(d);    

  t = zeros(ntriplets,1);                            
  for i = 1:ntriplets                       % Distances of central points along triplet paths
    t(i) = d(i)./sum(d(i:i+1));
  end;
  if (closed)
    t(ntriplets) = d(n-1)./sum(d([n-1,1]));
  end;

  xc = crds(:,1);
  yc = crds(:,2);

  xp = [xc(1:n-2) xc(2:n-1) xc(3:n)];       % Coordinates of triplets
  yp = [yc(1:n-2) yc(2:n-1) yc(3:n)];
  if (closed)
    xp = [xp; xc(n-1:n)' xc(2)];
    yp = [yp; yc(n-1:n)' yc(2)];
  end;

  xparam = zeros(ntriplets,3);              % Quadratic params for triplets
  yparam = zeros(ntriplets,3);
  for i = 1:ntriplets
    s = inv([1 0 0; 1 t(i) t(i)*t(i); 1 1 1]);
    xparam(i,:) = [s*xp(i,:)']';
    yparam(i,:) = [s*yp(i,:)']';
  end;

  nsp_orig = nsp;
  nsp = nsp*2;                              % Oversample the functions

  t1 = zeros(nsp,1);                        % First triplet on which point lies
  t2 = zeros(nsp,1);                        % Second triplet on which point lies
  d1 = zeros(nsp,1);                        % Distance along first triplet
  d2 = zeros(nsp,1);                        % Distance along second triplet

  if (closed)                               % Get relative dists for interpolated points
    [c,p] = pathpts(crds(1:n-1,:),nsp);     %   along path triplets
  else
    [c,p] = pathpts(crds,nsp);                 
  end;                                          

  % For each fitted point, find its relative positions 
  % along the two triplets on which it lies

  for i = 1:nsp-1
    j = max(find(cumpathlen<=p(i)));        % Find path node just before point
    if (j > 1)                              % If past second node,
      t2(i) = j-1;                          %   stash in first-triplet list
      d2(i) = (p(i)-cumpathlen(j-1)) / (cumpathlen(j+1)-cumpathlen(j-1));
    end;
    if (j==1 & closed)
      t2(i) = ntriplets;
      d2(i) = (p(i)+cumpathlen(ntriplets+1)-cumpathlen(ntriplets))...
              / (pathlen-cumpathlen(ntriplets));
    end;
    if (j <= ntriplets)
      t1(i) = j;                            % Stash in second-triplet list
      d1(i) = (p(i)-cumpathlen(j)) / (cumpathlen(j+2)-cumpathlen(j));           
    end;
  end;

  if (closed)
    t2(nsp) = ntriplets-1;
    d2(nsp) = 1;
    t1(nsp) = t2(1);
    d1(nsp) = d2(1);
  else
    t2(nsp) = ntriplets;
    d2(nsp) = 1;
  end;

  x = zeros(nsp,1);
  y = zeros(nsp,1);

  for i = 1:nsp                             % Get curve coords for each new point
    if (t2(i) == 0)                           % Point located in first internode 
      s = [1 d1(i) d1(i)*d1(i)]';             %   of open path
      x(i) = xparam(t1(i),:)*s;
      y(i) = yparam(t1(i),:)*s;

    elseif (t1(i) == 0)                       % Point located in last internode 
      s = [1 d2(i) d2(i)*d2(i)]';             %   of open path
      x(i) = xparam(t2(i),:)*s;
      y(i) = yparam(t2(i),:)*s;

    else                                      % Point in midst of open path or anywhere on closed path
      s = [1 d1(i) d1(i)*d1(i)]';             % Curve coords predicted by 1st part of 2nd triplet
      x1 = xparam(t1(i),:)*s;
      y1 = yparam(t1(i),:)*s;
      s = [1 d2(i) d2(i)*d2(i)]';             % Curve coords predicted by 2nd part of 1st triplet
      x2 = xparam(t2(i),:)*s;
      y2 = yparam(t2(i),:)*s;

      j = max(find(cumpathlen<=p(i)));        % Find path node just before point
      w =  (p(i)-cumpathlen(j)) / (cumpathlen(j+1)-cumpathlen(j));    % Weight

      x(i) = w*x1 + (1-w)*x2;
      y(i) = w*y1 + (1-w)*y2;
    end;
  end;

  x0 = crds(1,1);
  y0 = crds(1,2);
  if (closed)                                 % Extend curve to endpoints
    x = [x0; x; x0];
    y = [y0; y; y0];
  else
    xf = crds(n,1);
    yf = crds(n,2);
    x = [x0; x; xf];
    y = [y0; y; yf];
  end;

  scrds = [x y];                              % Concatenate into single matrix

  ns = ceil(0.03*length(x));            
  d1 = eucl(scrds(1:ns,:),scrds(1,:));        % Check behavior around endpoints
  dd = d1(2:ns) - d1(1:ns-1);
  i = max(find(dd<0));
  if (~isempty(i))
    scrds(2:i,:) = [];
  end;

  n = size(scrds,1);
  d2 = eucl(scrds([n:-1:n-ns+1],:),scrds(n,:));
  dd = d2(2:ns) - d2(1:ns-1);
  i = max(find(dd<0));
  if (~isempty(i))
    scrds([n-1:-1:n-i+1],:) = [];
  end;

  scrds = pathpts(scrds,nsp_orig);            % Resample the required number of pts

  if (doplot)                                 % Plot optional figure
    figure;
    plot(crds(:,1),crds(:,2),'ok');
    hold on;
    plot(scrds(:,1),scrds(:,2),'k');
    hold off;
    putbnd(scrds);
    axis('equal');
    axis('off');
  end;


  return;
