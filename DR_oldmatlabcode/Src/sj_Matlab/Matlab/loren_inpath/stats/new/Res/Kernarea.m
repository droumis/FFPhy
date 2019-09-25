% KERNAREA: Finds the area within a (1-alpha)-percent confidence region, 
%           based on the kernal density function of a 2D scatter of points.
%
%     Usage: [area,hglobal,h] = kernarea(crds,{alpha},{hglobal},{adapt},{doplot})
%
%         crds =    [n x 2] matrix of point coordinates.
%         alpha =   percentage for 1-alpha confidence region [default = 5].
%         hglobal = optional value of smoothing parameter; if not provided, a 
%                     maximum likelihood estimate is sought.
%         adapt =   optional boolean flag indicating, if true, that smoothing 
%                     parameters are to be 'adaptive', varying inversely in 
%                     response to local point density [default = 0].
%         doplot =  optional boolean flag indicating that a contour plot of the 
%                     probability density function is to be produced 
%                     [default = 0].
%         ---------------------------------------------------------------------
%         area =    area within the 1-alpha % confidence region.
%         hglobal = global smoothing parameter.
%         h =       [n x 1] vector of smoothing parameters for each observed point.
%

% RE Strauss, 5/11/00

function [area,hglobal,h] = kernarea(crds,alpha,hglobal,adapt,doplot)
  if (nargin < 2) alpha = []; end;
  if (nargin < 3) hglobal = []; end;
  if (nargin < 4) adapt = []; end;
  if (nargin < 5) doplot = []; end;

  if (isempty(alpha))
    alpha = 0.05;
  end;
  if (alpha > 1)
    alpha = alpha/100;
  end;

  [n,p] = size(crds);
  if (p ~= 2)
    error('  KERNAREA: 2D points only.');    
  end;

  xmean = mean(crds(:,1));              % Find bounds for square evaluation domain
  xmin = min(crds(:,1));
  xmax = max(crds(:,1));
  ymean = mean(crds(:,2));
  ymin = min(crds(:,2));
  ymax = max(crds(:,2));
  xrng = xmax - xmin;
  yrng = ymax - ymin;

  bounds = [xmin ymin; xmean ymin; xmax ymin; xmax ymean; 
            xmax ymax; xmean ymax; xmin ymax; xmin ymean];
  fx = ones(size(bounds,1),1);
  delta = 0.05;

  while (any(fx > alpha))               % Enlarge bounds till capture alpha contour
    [fx,hglobal,h] = kernpdf(bounds,crds,hglobal,adapt);

    if (any(fx([1 7 8]) > alpha))
      xmin = xmin - delta*xrng;
    end;
    if (any(fx([3 4 5]) > alpha))
      xmax = xmax + delta*xrng;
    end;
    if (any(fx([1 2 3]) > alpha))
      ymin = ymin - delta*yrng;
    end;
    if (any(fx([5 6 7]) > alpha))
      ymax = ymax + delta*yrng;
    end;

    xrng = xmax - xmin;
    yrng = ymax - ymin;
    bounds = [xmin ymin; xmean ymin; xmax ymin; xmax ymean; 
              xmax ymax; xmean ymax; xmin ymax; xmin ymean];
  end;

  npt = 30;
  x = linspace(xmin,xmax,npt);          % Create grid of evaluation pts
  y = linspace(ymin,ymax,npt);
  [X,Y] = meshgrid(x,y);
  gridpts = [X(:) Y(:)];

  fx = kernpdf(gridpts,crds,hglobal,adapt); % Evaluate density at grid pts
  FX = reshape(fx,npt,npt);

  cs = contourc(x,y,FX,[alpha alpha])';     % Find alpha contour(s)
  lencs = size(cs,1);
  area = 0;
  i = 1;
  while (i < lencs)                         % Cycle thru alpha contours
    npt = cs(i,2);
    perim = cs(i+1:i+npt,:);

    add_area = 0;
    j = 0;
    while(~add_area & j<n)                  % If contour contains at least 1 pt,
      j = j+1;                              %   add its area to total 
      add_area = isinpoly(crds(j,:),perim); %   else it is a valley rather than peak
    end;

    if (add_area)
      area = area + polyarea(perim);
    else
      area = area - polyarea(perim);
    end;

    i = i + npt + 1;
  end;

  if (doplot)
    [cs,ch] = contour(X,Y,FX,[0.05:0.05:max(fx)]);
    hold on;
    sqplot(gridpts);
    clabel(cs,ch);
    plot(crds(:,1),crds(:,2),'k*');
    hold off;
  end;

  return;

