% HULLPEEL: Finds a set of nested convex hulls characterizing 2-dimensional 
%           percentiles by hull peeling (Barnett 1976, Bebbington 1978).  
%           The (1-a)% hull is the largest hull containing no more 
%           than (1-a)% of the data (Riani Zani & Corbellini 1997).
%
%     Usage: [pt_pctile,hull_pctile,hull_area,med] = hullpeel(crds,{doplot})
%
%           crds =        [n x 2] matrix of point coordinates.
%           doplot =      optional boolean flag indicating, if true, that points 
%                           and hulls are to be plotted [default = 0].
%           --------------------------------------------------------------------
%           pt_pctile =   [n x 1] vector of percentiles of individual points.
%           hull_pctile = [h x 1] vector of percentile values of nested hulls,
%                           for h hulls.
%           hull_area =   [h x 1] vector of areas within each hull.
%           med =         [1 x 2] vector of coordinates of 2D median.
%

% Barnett, V. 1976. The ordering of multivariate data.  J.Am.Statist.Assoc. 
%   A139:318-339.
% Bebbington, AC. 1978. A method of bivariate trimming for robust estimation
%   of the correlation coefficient.  Appl.Statist. 27:221-226.
% Riani, M., S. Zani & A. Corbellini. 1997. Robust bivariate boxplots and
%   visualization of multivariate data.  In: I. Balderjahn, R. Mathar &
%   M. Schader (eds.), Classification, Data Analysis, and Data Highways,
%   pp. 93-100.  Springer-Verlag.

% RE Strauss, 1/5/99
%   5/31/99 - output the 2D median.
%   9/3/99 -  changed plot color for Matlab v5.
%   5/20/00 - added percentile values to plot; added hull areas.

function [pt_pctile,hull_pctile,hull_area,med] = hullpeel(crds,doplot)
  if (nargin<2) doplot = []; end;

  get_areas = 0;
  get_med = 0;
  if (nargout >= 3)
    get_areas = 1;
  end;
  if (nargout >= 4)
    get_med = 1;
  end;

  if (isempty(doplot))
    doplot = 0;
  end;
  if (doplot)
    get_med = 1;
  end;

  [n,p] = size(crds);
  if (p~=2)
    error('HULLPEEL: coordinates must be two-dimensional');
  end;

  hull_id = zeros(n,1);               % Allocate output matrices
  hull_pctile = [];
  hull_area = [];
  med = []; 

  if (doplot)                         % Initialize plot
    scatter(crds);
    hold on;
  end;

  incl = ones(n,1);
  h = 0;
  
  while (sum(incl)>=3)                 % Peel hulls
    p = sum(incl)*100/n;
    hull_pctile = [hull_pctile; p];
    in = find(incl==1);
    [hull_pts,index] = hull(crds(in,:));
    index = in(index);
    nh = length(index);
    incl(index) = zeros(nh,1);
    h = h+1;
    hull_id(index) = h*ones(nh,1);

    if (get_areas)
      hull_area = [hull_area; polyarea(hull_pts)];
    end;

    if (get_med)
      med = mean(crds(index,:));
    end;

    if (doplot)
      plot(hull_pts(:,1),hull_pts(:,2),'k');
      [mh,ih] = max(hull_pts(:,1));
      text(hull_pts(ih,1),hull_pts(ih,2),sprintf('  %1.0f',p));
    end;
  end;

  pt_pctile = hull_id;                % Identify percentiles of individual points
  len_hpctile = length(hull_pctile);  
  for h = 1:len_hpctile
    i = find(pt_pctile==h);
    pt_pctile(i) = hull_pctile(h)*ones(length(i),1);
  end;

  i = find(pt_pctile==0);             % Points within central hull
  if (~isempty(i))
    leni = length(i);
    p = leni*100/n;
    pt_pctile(i) = p*ones(length(i),1);

    if (get_med)
      if (leni>1)
        med = mean(crds(i,:));
      else
        med = crds(i,:);
      end;
    end;
  end;

  if (doplot)
    plot(med(1),med(2),'+k');
    hold off;
  end;

  return;