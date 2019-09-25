% GRIDVARY: Given a set of 2D point coordinates, finds the fit of quadrat counts 
%           to a negative-binomial, Poisson, or binomial distribution, 
%           varying the grid number.
%
%     Usage: [ngrid,X2,df,pr,distrib] = gridvary(pts,distrib,{doplot})
%
%           pts =     [n x 2] matrix of point coordinates.
%           distrib = vector of one or more distributions to be fitted:
%                       1 = negative binomial (mean < var)
%                       2 = poisson (mean = var)
%                       3 = binomial (mean > var)
%                           [default = all].
%           doplot =  optional boolean flag indicating, if true, that plots of 
%                       ngrid vs X2 and pr are to be produced [default = 0].
%           -------------------------------------------------------------------
%           ngrid =   vector of grid densities (quadrats per grid side).
%           X2 =      vector of chi-square statistic values.
%           df =      vector of corresponding degrees of freedom.
%           pr =      vector of corresponding chi-square probabilities.
%           distrib = distribution-id vector indicating distribution with which 
%                       preceding vector elements are associated.
%

% RE Strauss, 7/1/99
%   9/3/99 - plot changes for Matlab v5.

function [ngrid,X2,df,pr,distrib] = gridvary(pts,distrib,doplot)
  if (nargin < 2) distrib = []; end;
  if (nargin < 3) doplot = []; end;

  if (isempty(distrib))
    distrib = [1 2 3];
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  [n,p] = size(pts);
  if (p~=2)
    error('GRIDVARY: 2-dimensional coordinates only');
  end;

  if (doplot)
    scatter(pts);
    axis('equal');
  end;

  mincrd = min(pts)-eps;
  maxcrd = max(pts)+eps;

  ngrid = [];
  X2 = [];
  df = [];
  pr = [];
  dd = [];

  for di = 1:length(distrib)
    dist = distrib(di);

    ng = 1;
    maxcount = 2;

    while (maxcount > 1)
      ng = ng+1;                      % Increment quadrats/side
      counts = zeros(ng,ng);

      xlim = linspace(mincrd(1),maxcrd(1),ng+1);  % Determine lower quadrate bounds
      xlim(ng+1) = [];
      ylim = linspace(mincrd(2),maxcrd(2),ng+1);
      ylim(ng+1) = [];

      for i = 1:n                     % Add each point to a quadrate
        nx = max(find(pts(i,1)>xlim));
        ny = max(find(pts(i,2)>ylim));
        counts(nx,ny) = counts(nx,ny) + 1;
      end;
      maxcount = max(max(counts));

      if (maxcount > 1)
        if (dist==1)
          [k,mu,oP,eP,X2val,dfval,prval] = negbino(counts);   % Fit neg-binomial
          titl = 'Negative binomial';
        elseif (dist==2)
          [lambda,oP,eP,X2val,dfval,prval] = poisson(counts); % Fit Poisson
          titl = 'Poisson';
        elseif (dist==3)

          titl = 'Binomial';
        else
          error('GRIDVARY: invalid distribution identifier');
        end;

        ngrid = [ngrid; ng];
        X2 = [X2; X2val];
        df = [df; dfval];
        pr = [pr; prval];
        dd = [dd; dist];
      end;
    end;

    if (doplot)
      i = find(dd==dist);
      X2df = X2(i)./df(i);
      ngd = ngrid(i);
      prob = pr(i);

      figure;
      plot(ngd,X2df,'k');
      putbnd(ngd,X2df);
      puttick(ngd);
      putxlab('Number of quadrats / grid side');
      putylab('Chi-square value / df');
      puttitle(titl);

      figure;
      plot(ngd,prob,'k');
      hold on;
      plot([min(ngd) max(ngd)],[0.05 0.05],'k:');
      hold off;
      putbnd(ngd,prob);
      v = axis;
      axis([v(1:2) -0.05 1.05]);
      puttick(ngd);
      putxlab('Number of quadrats/grid side');
      putylab('Chi-square probability');
      puttitle(titl);
    end;
  end;

  distrib = dd;

  return;
