% RADIUSFN: 2D radius function for a closed counterclockwise polygon, as a 
%           function of length along the boundary (arc-length).  Optionally 
%           resamples the boundary by equal angular deviations about the center 
%           point.  Closes the polygon if it is open.
%
%     Usage: [r,theta,perim] = radiusfn(crds,{start},{center},{eqang},{noplot})
%
%           crds =    [n x 2] matrix of point coordinates.
%           start =   optional index (subscript) of starting point 
%                       [default = first point in 'crds'].
%           center =  optional center point [default = centroid].
%           eqang =   optional boolean flag indicating, if true, that the 
%                       polygon is to be resampled by equal angular deviations 
%                       about the center point [default = 0].
%           noplot =  optional boolean flag indicating, if true, that a plot of 
%                       the radius function is not to be produced [default = 0].
%           --------------------------------------------------------------------
%           r =       [n x 1] vector of radii.
%           theta =   angular deviations of radii from starting vector.
%           perim =   perimeter distances (arc-lengths) of radii from starting 
%                       point.
%

% RE Strauss, 3/5/00
%   3/21/00 - optionally resample by equal angular deviations.
%   3/12/01 - change plot option to noplot from doplot.

function [r,theta,perim] = radiusfn(crds,start,center,eqang,noplot)
  if (nargin < 2) start = []; end;
  if (nargin < 3) center = []; end;
  if (nargin < 4) eqang = []; end;
  if (nargin < 5) noplot = []; end;

  if (isempty(eqang))
    eqang = 0;
  end;
  if (isempty(noplot))
    noplot = 0;
  end;

  [n,p] = size(crds);
  if (p~=2)
    error('  RADIUSFN: 2-dimensional input coordinates only.');
  end;

  if (isempty(start))
    start = 1;
  end;
  if (~isscalar(start))
    error('  RADIUSFN: starting index must be a scalar.');
  elseif (start<1 | start>n)
    error('  RADIUSFN: starting index out of range.');
  end;

  if (isempty(center))                    % Center point
    [a,p,center] = polyarea(crds);
  else
    if (length(center)==1)                  % Index to point in 'crds'
      centerpt = crds(center,:);
      crds(center,:) = [];
      n = n-1;
    else
      if (length(center)~=2)
        error('  RADIUSFN: invalid center-point specification');
      end;
    end;
  end;

  if (start>1)                          % Shift points if necessary
    crds = [crds(start:n,:); crds(1:start-1,:)];
  end;
  if (eucl(crds(1,:),crds(n,:)) > (1e-6)*max(range(crds)))  % Close polygon
    crds = [crds; crds(1,:)];           % Append first point to end of sequence
    n = n+1;
  end;

  if (eqang)                            % Resample boundary by equal angles
    [crds,nonmono] = polyangl(crds,center);

    if (nonmono)
      disp('  RADIUSFN warning: Boundary not unique with respect to theta.');
      disp('                    Outer boundary intersections returned.');
    end;
  end;

  r = eucl(crds,center);                % Radii
  perim = zeros(n,1);

  Q = ones(n,1)*crds(1,:);
  P = ones(n,1)*center;
  theta = angl(Q,P,crds,1);             % Angles from first point
  theta(n) = 2*pi;

  for i = 2:n
    perim(i) = perim(i-1) + eucl(crds(i-1,:),crds(i,:));  % Accum perimeter length
  end;

  if (~noplot)
    axis('equal');
    hold on;
    plot(crds(1,1),crds(1,2),'k*');
    for i = 2:n
      plot(crds(i-1:i,1),crds(i-1:i,2),'k');
      plot([center(1),crds(i,1)],[center(2),crds(i,2)],'k');
    end;
    hold off;

    figure;
    subplot(2,1,1);
    plot(theta,r,'k');
    putbnd(theta,r);
    putxlab('Theta');
    putylab('Radius');

    subplot(2,1,2);
    plot(perim,r,'k');
    putbnd(perim,r);
    putxlab('Distance along boundary');
    putylab('Radius');
  end;

  return;

