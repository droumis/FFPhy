% QUADSPLINE: Fit a smoothed polynomial boundary to an open or closed 
%             (polygonal) 2-dimensional path.  The boundary is fitted using 
%             splined quadratic functions fitted to overlaping triplets 
%             of points, which are then smoothed linearly along the curve
%             (Overhauser curves).  The points returned are equally spaced along 
%             the boundary.  Does not use an interior center point.
%
%     Usage: scrds = quadspline(crds,{closed},{start},{noplot},{nsp})
%
%           crds =    [n x 2] set of point coordinates describing the path or 
%                       boundary.  If the first point is repeated as the last  
%                       point, the boundary is assumed to be closed.
%           closed =  optional boolean value indicating, if true, that
%                       boundary is to be treated as closed [default = 0];
%                       the boundary is also treated as closed if the first 
%                       point is repeated as the last.
%           start =   optional index (subscript) of starting point 
%                       [default = first point in 'crds'].
%           noplot =  optional boolean value indicating, if true, that plot of  
%                       points and fitted boundary is to be suppressed 
%                       [default=0].
%           nsp =     optional number of fitted points to be returned 
%                       [default=200].
%           ------------------------------------------------------------------
%           scrds =   [nsp x 2] matrix of smoothed fitted points.
%

% RE Strauss, 8/28 1998.
%   6/13/99 - add optional plot.
%   8/21/99 - changed plot colors for Matlab v5.
%   3/18/00 - changed output to single matrix rather than two vectors;
%             allowed for specification of start point.
%   3/26/00 - corrected overshot problem in neighborhoods of endpoints.
%   4/1/02 -  change name from quadspln();
%             change default to produce plot rather than suppress it;
%             add 'closebound' option;
%             do single check for closed boundary and set flag.

% Bowyer, A & J Woodwark. 1983. A programmer's geometry.  Butterworths, London.

function scrds = quadspline(crds,closed,start,noplot,nsp)
  if (nargin < 2) closed = []; end;
  if (nargin < 3) start = []; end;
  if (nargin < 4) noplot = []; end;
  if (nargin < 5) nsp = []; end;

  if (isempty(nsp))
    nsp = 200;
  end;
  if (isempty(noplot))
    noplot = 0;
  end;
  if (isempty(closed))
    closed = 0;
  end;
  if (isempty(start))
    start = 1;
  end;
  
  [n,p] = size(crds);
  distends = eucl(crds([1,n],:));           % Distance between endpoints
  distendtol = 0.01*max(range(crds));       % Tolerance for recognizing closed boundary
  if (distends < distendtol)
    closed = 1;
    crds = crds(1:n,:);
    n = n-1;
  end;
  firstpt = crds(1,:);

  if (p ~= 2)
    error('  QUADSPLINE: 2-dimensional coordinates only.');
  end;
  if (n < 4)
    error('  QUADSPLINE: minimum number of points is 4.');
  end;

  if (start<1 | start>n)
    error('  QUADSPLINE: start point out of range.');
  end;
  
  if (start > 1)                            % Shift around start point
    crds = [crds(start:n,:); crds(1:start-1,:)];  % Shift
  end;
  
  if (closed)
    crds = [crds([n-1,n,1:n,1:3],:)];
    n = size(crds,1);
  end;

  ntriplets = n-2;                          % Number of overlapping triplets of points

  d = crds(1:n-1,:) - crds(2:n,:);
  d = sqrt(d(:,1).*d(:,1)+d(:,2).*d(:,2));  % Distances between consecutive points
  cumpathlen = [0; cumsum(d)];              % Cum distances to points
  pathlen = sum(d);    

  t = zeros(ntriplets,1);                            
  for i = 1:ntriplets                       % Distances of central points along triplet paths
    t(i) = d(i)./sum(d(i:i+1));
  end;

  xc = crds(:,1);
  yc = crds(:,2);
  xp = [xc(1:n-2) xc(2:n-1) xc(3:n)];       % Coordinates of triplets
  yp = [yc(1:n-2) yc(2:n-1) yc(3:n)];

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

  [c,p] = pathpts(crds,nsp);                % Get relative dists for interpolated points
                                            %   along path triplets

  % For each fitted point, find its relative positions 
  % along the two triplets on which it lies

  for i = 1:nsp-1
    j = max(find(cumpathlen<=p(i)));        % Find path node just before point
    if (j > 1)                              % If past second node,
      t2(i) = j-1;                          %   stash in first-triplet list
      d2(i) = (p(i)-cumpathlen(j-1)) / (cumpathlen(j+1)-cumpathlen(j-1));
    end;
    if (j <= ntriplets)
      t1(i) = j;                            % Stash in second-triplet list
      d1(i) = (p(i)-cumpathlen(j)) / (cumpathlen(j+2)-cumpathlen(j));           
    end;
  end;

  t2(nsp) = ntriplets;
  d2(nsp) = 1;

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
  
  if (~closed)
    x = [crds(1,1); x; crds(n,1)];            % Extend curve to endpoints
    y = [crds(1,2); y; crds(n,2)];
  end;
  scrds = [x y];                              % Concatenate into single matrix

  if (closed)                                 % Trim excessive ends for closed form
    d = eucl(scrds,firstpt);
    incl = ones(size(d));
    dprev = 1e30;
   
    i = 1;
    while (d(i)<=dprev)
      incl(i) = 0;
      dprev = d(i);
      i = i+1;
    end; 
    
    dprev = 1e30;
    n = size(scrds,1);
    i = n;
    while (d(i)<=dprev)
      incl(i) = 0;
      dprev = d(i);
      i = i-1;
    end;
    
    scrds = [firstpt; scrds(incl>0,:); firstpt];
  end;

  scrds = pathpts(scrds,nsp_orig);            % Resample the required number of pts

  if (~noplot)                                % Plot figure
    figure;
    plot(crds(3:end-2,1),crds(3:end-2,2),'ok',...
         crds(3,1),crds(3,2),'*k',...
         scrds(:,1),scrds(:,2),'k');
    putbnd(scrds);
    axis('equal');
    axis('off');
  end;

  return;
