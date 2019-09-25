% KERNPDF: Estimates adaptive 2D kernal density estimates.
%
%     Usage: [fx,h,hglobal] = kernpdf(x,xi,{hglobal},{adapt},{doplot})
%
%         x =       [m x 2] matrix of coordinates for points at which function 
%                     is to be evaluated.
%         xi =      [n x 2] matrix of observed point coordinates.
%         hglobal = optional value of smoothing parameter; if not provided, a 
%                     maximum likelihood estimate is sought.
%         adapt =   optional boolean flag indicating, if true, that smoothing 
%                     parameters are to be 'adaptive', varying inversely in 
%                     response to local point density [default = 0].
%         doplot =  optional boolean flag indicating that a contour plot of the 
%                     probability density function is to be produced 
%                     [default = 0].
%         ---------------------------------------------------------------------
%         fx =      [m x 1] vector of density estimates corresponding to x.
%         h =       [n x 1] vector of smoothing parameters for each observed point.
%         hglobal = global smoothing parameter.
%

% Worton, BJ. 1989. Kernel methods for estimating the utilization distribution 
%   in home-range studies. Ecology 70:164-168.
% Brunsdon, C. 1995. Estimating probability surfaces for geographic point data: 
%   an adaptive kernel algorithm. Computers & Geosciences 21:877-894.

% RE Strauss, 5/6/00

function [fx,h,hglobal] = kernpdf(x,xi,hglobal,adapt,doplot)
  if (nargin < 3) hglobal = []; end;
  if (nargin < 4) adapt = []; end;
  if (nargin < 5) doplot = []; end;

  if (isempty(adapt))
    adapt = 0;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  [n,pn] = size(xi);
  [m,pm] = size(x);
  if (pn~=2 | pm~=2)
    error('  KERNDADAPT: observed and evaluation points must be 2D.');
  end;

  if (isempty(hglobal))
    d = eucl(xi);
    d = trilow(d);
    hmin = min(d);
    hmax = mean(d);

    hglobal = fmin('kernlh',hmin,hmax,[],xi); % Find max-likelihood h

    if (doplot)
      hval = linspace(hmin,hmax);
      lk = zeros(size(hval));
      for i = 1:length(hval)
        lk(i) = -kernlh(hval(i),xi);
      end;
      figure;
      plot(hval,lk);
      putxlab('h');
      putylab('Log-likelihood');
    end;
  end;

  if (adapt)
    fx = kernden(xi,xi,hglobal);              % Evaluate observed points
    g = exp(mean(log(fx)));
    h = hglobal * ((fx./g).^(-0.5));           % h's at observed points 
  else
    h = hglobal;
  end;

  fx = kernden(x,xi,h);                     % Densities at test points

  if (doplot)
    n = 40;
    r = range([x;xi])/2;
    xmin = min([x;xi]) - r;
    xmax = max([x;xi]) + r;

    [X,Y] = meshgrid(linspace(xmin(1),xmax(1),n),linspace(xmin(2),xmax(2),n));
    z = kernden([X(:) Y(:)],xi,h);
    Z = reshape(z,n,n);

    figure;
    [cs,ch] = contour(X,Y,Z);
    clabel(cs,ch);
    axis('square');
    axis([xmin(1) xmax(1) xmin(2) xmax(2)]);
    hold on;
    plot(xi(:,1),xi(:,2),'k*');
    hold off;
  end;

  return;


