% CUBESPLN: Cubic-spline smoothing of a 2D boundary specified by a series of 
%           point coordinates.  Uses the radius as a function of angle, 
%           so an interior center point must be specified (defaults to centroid) 
%           and the boundary function must be unique with respect to angle.
%
%     Usage: [scrds,r,theta] = cubespln(crds,{start},{center},{doplot},{nsp})
%
%           crds =      [n x 2] matrix of point coordinates.
%           start =     optional index (subscript) of starting point 
%                         [default = first point in 'crds'].
%           center =    optional center point [default = centroid].
%           doplot =    optional boolean flag indicating, if true, that a plot 
%                         of the coordinates and spline is to be produced 
%                         [default = 0].
%           nsp =       optional number of points describing the spline 
%                         [default = 200].
%           ------------------------------------------------------------------
%           scrds =     [nsp x 2] set of cartesian spline coordinates.
%           [r,theta] = vectors (length nsp) of polar spline coordinates.
%

% RE Strauss, 3/5/00
%   3/19/00 - add optional start and center points.
%   9/12/01 - fixed bug with radius function.

function [scrds,r,theta] = cubespln(crds,start,center,doplot,nsp)
  if (nargin < 2) start = []; end;
  if (nargin < 3) center = []; end;
  if (nargin < 4) doplot = []; end;
  if (nargin < 5) nsp = []; end;

  if (isempty(nsp))
    nsp = 200;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;
  if (isempty(start))
    start = 1;
  end;

  [n,p] = size(crds);
  if (p~=2)
    error('  CUBESPLN: 2-dimensional input coordinates only.');
  end;

  if (isempty(center))                    % If center point not given,
    [a,p,center] = polyarea(crds);        %   get centroid
  else
    center = crds(center,:);              % Else save point crds
    crds(center,:) = [];                  %   and remove point from list
    n = n-1;
  end;

  if (start<1 | start>n)
    error('  CUBESPLN: start point out of range.');
  end;

  if (eucl(crds([1,n],:)) < (1e-6)*max(range(crds)))  % Open the polygon
    crds(n,:) = [];
    n = n-1;
  end;

  if (start > 1)                          % Shift around start point
    crds = [crds(start:n,:); crds(1:start-1,:)];  % Shift
  end;

  crds = [crds; crds; crds];              % Replicate coordinates

  theta0 = angledev(center+[-1 0],center,crds(1,:));  % Ref angle for 1st crd

  [r,theta] = radiusfn(crds,[],[],[],1);  % Radius function
  lentheta = length(theta);
  d = theta(2:lentheta) - theta(1:lentheta-1);
  while (any(d<0))
    i = find(d<0);
    i = i(1)+1;
    theta(i:lentheta) = theta(i:lentheta) + 2*pi;
    d = theta(2:lentheta) - theta(1:lentheta-1);
  end;

  ti = linspace(theta(1),theta(length(theta)),3*nsp); % Cubic spline interpolation
  ri = interp1(theta,r,ti,'spline');

  i = min(find(ti>=theta(1)+2*pi));       % Pull out middle third of interp pts
  j = max(find(ti<=theta(1)+4*pi));

  theta = linspace(ti(i),ti(j+1),nsp)';   % Reinterpolate
  r = interp1(ti',ri',theta);
  theta = theta + theta0;

  [x,y] = polarcrd(r,theta,1);            % Retransform to cartesian crds
  scrds = [x y] + ones(nsp,1)*center;

  if (doplot)
    figure;
    plot(crds(:,1),crds(:,2),'ko');
    axis('square');
    axis(sqplot(crds));
    hold on;
    plot(scrds(:,1),scrds(:,2),'k');
    hold off;
  end;

  return;

