% HOMERANGE:  Estimates home range for various measures, based on spatial 
%             coordinates of points of capture.
%
%     Usage: [hr,hrci] = homerange(crds,{stat},{iter},{ci_level})
%
%             crds =      [n x 2] matrix of coordinates of points of capture.
%             stat =      statistic used to estimate home range:
%                           1 = area of convex hull [default];
%                           2 = mean squared deviation from centroid;
%                           3 = mean squared deviation from 50% trimmed centroid;
%                           4 = median squared deviation from medioid.
%             iter =      optional number of bootstrap iterations [default = 0].
%             ci_level =  optional level of bootstrapped confidence interval 
%                           [default = 95].
%             -------------------------------------------------------------------
%             hr =        home-range estimated (squared units).
%             hrci =      upper and lower bounds [lb,ub] of confidence interval.
%                       

% RE Strauss, 7/26/00

function [hr,hrci] = homerange(crds,stat,iter,ci_level)
  if (nargin < 2) stat = []; end;
  if (nargin < 3) iter = []; end;
  if (nargin < 4) ci_level = []; end;

  [n,p] = size(crds);
  if (p~=2)
    error('  HOMERANGE: two-dimensional coordinates only.');
  end;

  if (isempty(stat))
    stat = 1;
  end;
  if (isempty(iter))
    iter = 0;
  end;
  if (isempty(ci_level))
    ci = 0.95;
  end;

  hr = homerangef(crds,[],[],[],stat);

  hrci = [];
  if (iter)
    hrci = bootstrp('homerangef',1,iter,1-ci_level,crds,[],[],stat);
  end;

  return;
