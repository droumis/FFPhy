% POLYANGL: Resamples the boundary of a polygon by equal angles about a central 
%           point (the centroid, by default).  If the boundary coordinates are 
%           given in a clockwise direction, the angles sampled are clockwise; 
%           and vice versa.  If the input polygon is closed, the first and last 
%           resampled points will be identical; if the polygon is open, the 
%           resampled polygon will also be open.  If the boundary points are 
%           not monotonic with respect to angle from center, a 'nonmono' flag 
%           is returned.  If a radius from center intersects the boundary in 
%           more than one place due to non-monotonicity, the outermost boundary 
%           intersection is returned.
%             Optionally allow set of radii angles (radians) to be passed for 
%           resampling.  The function then returns the points of intersection 
%           of the polygon boundary with these radii.
%
%     Usage: [acrds,nonmono] = polyangl(crds,{center},{angles},{npts})
%
%           crds =    [n x 2] matrix of point coordinates of boundary.
%           center =  optional 2-element vector of center coordinates 
%                       [default = centroid].
%           angles =  optional set of radii angles (assumed counterclockwise) 
%                       for resampling.
%           npts =    optional number of resampled boundary points to be 
%                       returned [default = n or length of 'angles'].
%           -----------------------------------------------------------------------
%           acrds =   [npts x 2] resampled boundary points
%           nonmono = boolean flag indicating, if true, that the boundary 
%                       points are not monotonic with respect to angle from center.
%           

% RE Strauss, 3/25/00
%   3/26/00 - added option to pass radii angles.

function [acrds,nonmono] = polyangl(crds,center,angles,npts)
  if (nargin < 2) center = []; end;
  if (nargin < 3) angles = []; end;
  if (nargin < 4) npts = []; end;

  [n,p] = size(crds);
  if (p ~= 2)
    error('  POLYANGL: 2-dimensional coordinates only.');
  end;

  if (isempty(npts))
    npts = n;
  end;

  angles_passed = 0;
  if (isempty(angles))
    rtheta = linspace(0,2*pi,npts)';        % Equal-dev angles
  else
    rtheta = angles;                        % Angles given as input
    npts = length(angles);
    angles_passed = 1;
  end;

  polyopen = 0;
  if (eucl(crds(1,:),crds(n,:)) > 1e-6)     % Close polygon
    crds = [crds; crds(1,:)];
    n = n+1;
    if (~angles_passed)
      npts = npts+1;
    end;
    polyopen = 1;
  end;

  if (isempty(center))
    [a,p,center] = polyarea(crds);
  end;

  if (length(center(:)) ~= 2)
    error('  POLYANGL: invalid center point.');
  end;

  savecrds = crds;
  reftheta = angl(center+[1 0],center,crds(1,:),1);  % Rotate config
  crds = rotate(crds,-reftheta,center);
  if (angles_passed)
    rtheta = rtheta - reftheta;
    i = find(rtheta<0);
    if (~isempty(i))
      rtheta(i) = rtheta(i) + 2*pi;
    end;
  end;
  
  crd_theta = angl(ones(n,1)*(center+[1 0]),center,crds,1); % Get angles for pts

  reflect = 0;                            % Reflect crds if not counterclockwise
  ncc = sum(crd_theta(2:n) > crd_theta(1:n-1));  
  if (ncc < n/2)
    crds = savecrds;                        % Restore
    crds(:,1) = -crds(:,1);                 % Reflect
    center(1) = -center(1);
    reftheta = angl(center+[1 0],center,crds(1,:),1);  % Rotate reflected config
    crds = rotate(crds,-reftheta,center);
    crd_theta = angl(ones(n,1)*(center+[1 0]),center,crds,1); % Get angles for pts
    reflect = 1;
  end;
  crd_theta([1 n]) = [0 2*pi];

  nonmono = 0;                            % Determine whether bounary pts are
  ncc = sum(crd_theta(2:n) < crd_theta(1:n-1)); % monotonic with respect to angle
  if (ncc)                              
    nonmono = 1;
  end;

  bound = [crds(1:n-1,:) crds(2:n,:)];    % Boundary line segments

  r = 1.5*max(range(crds));               % Length of radius
  if (angles_passed)
    [x,y] = polarcrd(r,rtheta,1);
    radii = [ones(npts,1)*center [x y] + ones(npts,1)*center];
  else
    ccrds = circcrds(r,center,npts);        % Endpoints of radii
    radii = [ones(npts,1)*center ccrds];    % Radii line segments
  end;

  acrds = [crds(1,:); zeros(npts-2,2); crds(1,:)];  % Allocate output matrix

  for ir = 2:(npts-1)                     % Get intersections of radii with boundary
    i = max(find(crd_theta <= rtheta(ir)));
    x = [];
    while (isempty(x))
      [intsct,x,y] = intrsect(bound(i,:),radii(ir,:));  % Get intersection
      i = i+1;
    end;
    acrds(ir,:) = [x y];
  end;

  acrds = rotate(acrds,reftheta,center);  % Rotate to original orientation
  if (reflect)                            % Re-reflect, if necessary
    acrds(:,1) = -acrds(:,1);                 
  end;
  if (polyopen)                           % Open polygon, if necessary
    acrds(npts,:) = [];
  end;

  return;
