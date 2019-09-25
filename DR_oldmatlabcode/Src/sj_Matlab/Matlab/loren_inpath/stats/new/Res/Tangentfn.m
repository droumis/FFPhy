% TANGENTFN: 2D tangent-angle function of a quadratic-spline smoothing of a 
%           closed polygon, as a function of length along the boundary 
%           (arc-length).  Closes the polygon if it is open.  Tangent angle of 
%           the start point is defined to be zero.  Optionally plots the fitted 
%           outline and the tangent-angle function.
%
%     Usage: [theta,perim] = tangentfn(crds,{start},{noplot},{nsp})
%
%           crds =    [n x 2] matrix of point coordinates.
%           start =   optional index (subscript) of starting point 
%                       [default = first point in 'crds'].
%           noplot =  optional boolean flag indicating, if true, that plot of 
%                       the radius function is not to be produced [default = 0].
%           nsp =     optional number of fitted points to be returned 
%                       [default=200].
%           -------------------------------------------------------------------
%           theta =   [n x 1] vector of tangent-angles.
%           perim =   perimeter distances (arc-lengths) of tangent-angles 
%                       from starting point.
%

% RE Strauss, 3/26/00
%   9/12/01 - added initial figure window.

function [theta,perim] = tangentfn(crds,start,noplot,nsp)
  if (nargin < 2) start = []; end;
  if (nargin < 3) noplot = []; end;
  if (nargin < 4) nsp = []; end;

  get_perim = 0;
  if (nargout > 1)
    get_perim = 1;
  end;

  if (isempty(start))
    start = 1;
  end;
  if (isempty(noplot))
    noplot = 0;
  end;
  if (isempty(nsp))
    nsp = 200;
  end;

  if (noplot)
    get_perim = 1;
  end;

  [n,p] = size(crds);
  if (p~=2)
    error('  TANGENTFN: 2-dimensional input coordinates only.');
  end;

  if (~isscalar(start))
    error('  TANGENTFN: starting index must be a scalar.');
  elseif (start<1 | start>n)
    error('  TANGENTFN: starting index out of range.');
  end;

  if (eucl(crds(1,:),crds(n,:)) > (1e-6)*max(range(crds)))  % Close polygon
    crds = [crds; crds(1,:)];           % Append first point to end of sequence
    n = n+1;
  end;

  scrds = quadspln(crds,start,0,nsp);   % Quadratic splining of crds
  scrds = [scrds(nsp-1,:); scrds; scrds(2,:)];  % Extend fitting

  P = scrds(1:nsp,:);                   % Points forming angles
  Q = P + ones(nsp,1)*([1 0]);
  R = scrds(3:nsp+2,:);

  theta = angl(Q,P,R,1);                % Tangent-angles
  theta(nsp) = theta(1) + 2*pi;
  dt = [theta(2:nsp) - theta(1:nsp-1)];

  i = find(dt > pi);
  if (~isempty(i))
    dt(i) = dt(i) - 2*pi;
  end;
  i = find(dt < -pi);
  if (~isempty(i))
    dt(i) = dt(i) + 2*pi;
  end;
  theta = [0; cumsum(dt)];

  if (get_perim)                        % Distance along perimeter
    perim = zeros(nsp,1);
    for i = 2:nsp
      perim(i) = perim(i-1) + eucl(scrds(i,:),scrds(i+1,:));  % Accum perimeter length
    end;
  end;

  if (~noplot)
    figure;
    plot(crds(:,1),crds(:,2),'ko');
    hold on;
    plot(crds(start,1),crds(start,2),'k*');
    plot(scrds(:,1),scrds(:,2),'k');
    hold off;
    putbnd([crds; scrds]);

    figure;
    plot(perim,theta,'k');
    putbnd(perim,theta);
    putxlab('Distance along boundary');
    putylab('Tangent-angle');
  end;

  return;


