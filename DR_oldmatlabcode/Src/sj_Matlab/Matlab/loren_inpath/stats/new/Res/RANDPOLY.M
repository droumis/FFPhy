% RANDPOLY: Generates uniform-random points within a given arbitrary polygon.  
%           Use RANDCNVX if the polygon is convex.
%
%     Syntax: pt = randpoly(poly,{npts})
%
%           poly - [m x 2] matrix of coordinates specifying a polygon.
%           npts - number of random points to be generated [default=1].
%           -----------------------------------------------------------
%           pt -   [npts x 2] matrix of random point coordinates.
%

% RE Strauss, 7/17/96
%   9/7/99 - changed plot colors for Matlab v5.

function pt = randpoly(poly,npts)
  if (nargin < 2) npts = []; end;

  if (isempty(npts))                  % Default number of points
    npts = 1;
  end;

  [m,p] = size(poly);
  if (p~=2)                           % Check dimensions
    error('  Procedure implemented only for 2 dimension');
  end;
  if (poly(m,:)~=poly(1,:))           % Close polygon if necessary
    poly = [poly; poly(1,:)];
    m = m+1;
  end;

  xmin = min(poly(:,1));              % Rectangle enclosing polygon
  xmax = max(poly(:,1));
  xrange = xmax - xmin;
  ymin = min(poly(:,2));
  ymax = max(poly(:,2));
  yrange = ymax - ymin;

  pt = zeros(npts,2);
  initpt = 1;
  n = npts;

  while (n>0)
    p = [rand(n,1)*xrange+xmin, rand(n,1)*yrange+ymin];

    isin = isinpoly(p,poly);
    indx = find(isin);

    finalpt = initpt + length(indx) - 1;
    pt(initpt:finalpt,:) = p(indx,:);

    initpt = finalpt+1;
    n = n-length(indx);
  end;

% clf;plot(poly(:,1),poly(:,2),'-',pt(:,1),pt(:,2),'ko');

  return;

