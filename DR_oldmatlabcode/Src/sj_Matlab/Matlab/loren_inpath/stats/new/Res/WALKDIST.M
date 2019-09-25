% WALKDIST: Models a random walk, in 2 dimensions, of two descendant taxa
%           from an ancestral taxon positioned at the origin of the coordinate
%           system.  Plots corresponding pairwise distances traveled, as well
%           as specified percentile regions.
%
%     Usage: [dist1,dist2] = walkdist(iter,nsteps,stepflag,perconf)
%
%           iter =     number of iterations of paired walks [default = 100].
%           nsteps =   number of steps      [default = 25].
%           stepflag = flag indicating distribution from which step-lengths
%                        are sampled:
%                         0 = normal,  abs(N(0,1)) [default]
%                         1 = uniform, U[0,1]
%                         2 = constant 1
%           perconf =  vector of percentile confidence values to be plotted
%                        as rays from the origin [default: 5,25,50,75,95].
%           ----------------------------------------------------------------
%           dist1 =    vector of Euclidean distances from origin of taxon 1
%                        at end of walks.
%           dist2 =    vector of Euclidean distances from origin of taxon 2
%                        at end of walks.
%

% RE Strauss, 8/13/95
%   8/21/99 - changed plot colors, and other misc statements, for Matlab v5.

function [dist1,dist2] = walkdist(iter,nsteps,stepflag,perconf)
  if (nargin < 1) iter = []; end;
  if (nargin < 2) nsteps = []; end;
  if (nargin < 3) stepflag = []; end;
  if (nargin < 4) perconf = []; end;

  if (isempty(iter))
    iter = 100;
  end;
  if (isempty(nsteps))
    nsteps = 25;
  end;
  if (isempty(stepflag))
    stepflag = 0;
  end;
  if (isempty(perconf))
    perconf = [5 25 50 75 95];
  end;

  dist1 = zeros(iter,1);              % Allocate output distance vectors
  dist2 = dist1;
  theta = dist1;                      % Allocate
  perconf = sort(perconf')';          % Sort percentiles into ascending seq

  for it = 1:iter
    dist1(it) = walkrand(2,nsteps,stepflag);  % Final distance for taxon 1
    dist2(it) = walkrand(2,nsteps,stepflag);  % Final distance for taxon 2
    theta(it) = dist2(it)/dist1(it);          % Tangent-angle of ray from origin
  end;

  [t,indx] = sort(theta);             % Sort tangent-angles
  m = ceil(max(max([dist1 dist2])));  % Find max value for plot

  clf;
  hold on;
  axis([0 m 0 m]);
  axis('square');

  plot(dist1,dist2,'ob');             % Plot pairwise distances
  plot([0 m],[0 m],'k-');             % Plot line of equal distances from anc
  text(m+0.015*m, m, 'D1=D2');

  for i = 1:length(perconf)           % Find and plot confidence rays
    pos = iter * perconf(i) / 100;
    a = theta(indx(pos));             % Find point at percentile ray
    x = m;                            % Coordinates of ray
    y = a*x;
    if (y>m)
      y = m;
      x = y/a;
    end;
    plot([0 x],[0 y],'r-');
    text(x+0.015*m, y, num2str(perconf(i)));
  end;

  xlabel('Distance from ancestor to T1');
  ylabel('Distance from ancestor to T2');
  hold off;

  return;

