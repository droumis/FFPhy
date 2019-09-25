% THINPLATE: Thin-plate splines and principal warps.
%
%     Usage: [warp,eval,evec,D,theta1,theta2] = ...
%                   thinplate(base,targ,{keepscale},{gridsize},{trimgrids})
%
%         base =      [n x 2] matrix of base points.
%         targ =      [n x 2] matrix of correponding target points.
%         keepscale = optional boolean flag indicating whether (=1) or not (=0) 
%                       to maintain original scaling of base and target point 
%                       configurations on grid displays [default=1].
%         gridsize =  optional number of points along edge of square grid for display
%                       [default=20].
%         trimgrids = optional boolean flag indicating that base and target grids are
%                       to be trimmed within the convex hulls of the landmarks
%                       [default = 0].
%         ---------------------------------------------------------------------------
%         warp =      [n x n] matrix of bending energy as a function of changes
%                       in the target coordinates.
%         eval =      [k x 1] vector of k<n nonzero eigenvalues of warp.
%         evec =      [n x k] matrix of corresponding eigenvectors (columns).
%         D =         principal dilatations of affine transformation.
%         theta1 =    directions (degrees counterclockwise) of principal dilatations.
%         theta2 =    final rotation (degrees counterclockwise) of affine 
%                       transformation.
%

% RE Strauss, 8/20/98, modified from Eric Dyreson's SAS/IML code, based on Bookstein 1989.
%   8/20/99 - changed plot colors for Matlab v5.
%   1/4/00 -  changed usage of sqplot().
%   4/6/02 -  change ellips() to ellipsebound();
%             add error message for forms with different numbers of landmarks;
%             rename from thinplat() to thinplate();
%             add 'trimgrids' option.

% Bookstein, F.L.  1989.  Principal warps: Thin-plate splines and the
%   decomposition of deformations.  IEEE Trans. Patt. Anal. Mach. Intell.
%   11:567-585.

function [warp,eval,evec,D,theta1,theta2] = ...
                  thinplate(base,targ,keepscale,gridsize,trimgrids)
  if (nargin < 3) keepscale = []; end;
  if (nargin < 4) gridsize = [];  end;
  if (nargin < 5) trimgrids = []; end;

  if (isempty(keepscale)) keepscale = 1; end;
  if (isempty(trimgrids)) trimgrids = 0; end;

  if (isempty(gridsize))  
    steps = 19;
  else
    steps = gridsize-1;
  end;

  n = size(base,1);
  if (size(targ,1) ~= n)
    error('  THINPLATE: base and target forms must have same number of landmarks.');
  end;

  base = zcenter(base);             % Center point configs on origin
  targ = zcenter(targ);

  x = (base(:,1)*ones(1,n) - ones(n,1)*base(:,1)').^2;
  y = (base(:,2)*ones(1,n) - ones(n,1)*base(:,2)').^2;
  r = sqrt(x+y);                      % U function evaluation
  r2 = r.^2;
  logr2 = log(r2 + eye(n));

  K = r2 .* logr2;                    % Construct L matrix
  P = [ones(n,1) base];
  L = [K P; P' zeros(3,3)];

  Linv = inv(L);                      % Warp matrix
  Lred = Linv(1:n,1:n);
  warp = Lred * K * Lred;

  LV = Linv * [targ; zeros(3,2)];     % Transformation map
  A = LV((n+2):(n+3),:)';             % Linear terms
  [U,S,V] = svd(A);                   % Decomposition
  D = diag(S);                        % Affine dilatations
  theta1 = zeros(2,1);                % Orientations of dilatations
  theta1(1) = -asin(V(1,2));
  theta1(2) = theta1(1)+pi/2;
  t = theta1(1);
  T1 = [cos(t) sin(t);-sin(t) cos(t)]; % Orientation matrix
  theta2 = asin(U(1,2))+t;             % Net affine rotation
  t = -theta2;
  T2 = [cos(t) sin(t);-sin(t) cos(t)]; % Net rotation matrix

  [evec,eval] = eigen(warp);          % Eigenanalysis of warp matrix
  i = find(eval>eps);
  eval = eval(i);
  evec = evec(:,i);

  bounds = sqplot(base,[],1);         % Make rectangular grid of points
  xmin = bounds(1);
  xmax = bounds(2);
  ymin = bounds(3);
  len_step = (xmax-xmin)/steps;
  m1 = meshgrid(0:steps,0:steps);
  m2 = m1';
  intgrid = [m1(:) m2(:)];
  sqgrid = [intgrid(:,1)*len_step+xmin, intgrid(:,2)*len_step+ymin];

  Tgrid = [];                         % Non-affine transformation mapping
  for i=1:size(sqgrid,1)
    U = ((base-ones(n,1)*sqgrid(i,:)).^2)*[1; 1];
    U = U .* log(U+(U==0));
    U = [U; 1; sqgrid(i,:)']';
    Tgrid = [Tgrid; U*LV];
  end;
  Tgrid = Tgrid * T2;                 % Affine rotation for display
  targ = targ * T2;

  [map1,map2] = lstra(base,targ);     % Orthogonal least-squares Procrustes mapping

  if (keepscale)
    xmin = min(min([sqgrid(:,1) Tgrid(:,1)]));
    xmax = max(max([sqgrid(:,1) Tgrid(:,1)]));
    ymin = min(min([sqgrid(:,2) Tgrid(:,2)]));
    ymax = max(max([sqgrid(:,2) Tgrid(:,2)]));
  end;
  
  if (trimgrids)
    isin = isinpoly(sqgrid,hull(base));
    sqgrid = sqgrid(isin==1,:);
    isin = isinpoly(Tgrid,hull(targ));
    Tgrid = Tgrid(isin==1,:);
  end;

  clf;
  subplot(2,2,1);
  hold on;
  title('Reference form');
  plot(sqgrid(:,1),sqgrid(:,2),'b.');
  plot(base(:,1),base(:,2),'r*');
  if (keepscale)
    axis([xmin xmax ymin ymax]);
  else
    axis([min(sqgrid(:,1)) max(sqgrid(:,1)) min(sqgrid(:,2)) max(sqgrid(:,2))]);
  end;
  axis('square');
  axis('equal');
  axis('off');
  hold off;

  subplot(2,2,2);
  hold on;
  title('Target form');
  plot(Tgrid(:,1),Tgrid(:,2),'b.');
  plot(targ(:,1),targ(:,2),'ro');
  if (keepscale)
    axis([xmin xmax ymin ymax]);
  else
    axis([min(Tgrid(:,1)) max(Tgrid(:,1)) min(Tgrid(:,2)) max(Tgrid(:,2))]);
  end;
  axis('square');
  axis('off');
  hold off;

  subplot(2,2,3);
  hold on;
  title('Affine dilatations');
  a = D(1)/2;                         % Half-dilatations
  b = D(2)/2;
  t = theta1(1);                      % Orientation of major dilatation
  [x,y] = ellipsebound(a,b,0,0,t);    % Locus of ellipse
  plot(x,y,'k');                      % Plot ellipse
  ma = [-a 0;a 0]*T1;                 % Plot major axis
  plot(ma(:,1),ma(:,2),'k');
  ma = [0 -b;0 b]*T1;                 % Plot minor axis
  plot(ma(:,1),ma(:,2),'k');
  bounds = sqplot(x,y,40);
  axis(bounds);
  axis('square');
  axis('off');
  hold off;

  subplot(2,2,4);
  hold on;
  title('Procrustes mapping');
  plot(map1(:,1),map1(:,2),'k*');     % Plot base coordinates
  plot(map2(:,1),map2(:,2),'ko');     % Plot superimposed target coordinates
  axis('equal');
  axis('square');
  axis('off');
  hold off;

  theta1(1) = theta1(1)*360/(2*pi);   % Convert angles to degrees
  theta1(2) = theta1(1)+90;
  theta2 = theta2*360/(2*pi);

  return;
