% MSTCLUST: Find 'significant' clusters of points in 2-space based on
%           deletion of a statistically significant subset of longest edges of 
%           the minimum-spanning tree, with respect to randomized distributions 
%           of points within the convex hull of the data points.
%
%     Syntax: [ngrps,grps,prob] = mstclust(crds,{iter},{critval},{doplot})
%
%           crds -    [n x 2] matrix of point coordinates.
%           iter -    number of randomization iterations [default=1000].
%           critval - critical tail-value for deleting longest edges 
%                       [default=0.15].
%           doplot -  boolean value indicating whether (1, =TRUE) 
%                       or not (0, =FALSE) scatter plots of points, hulls, and 
%                       trees are to be produced [default = 0].
%           ------------------------------------------------------------------
%           ngrps -   number of groups recognized at the critical value.
%           grps -    [n x 1] vector of assigned group membership for all points.
%           prob -    [(n-1) x 1] vector of randomized tail probabilities for 
%                       observed M-S tree edges (branches), ranked longest to shortest.
%

% RE Strauss, 12/31/96
%   9/7/99 - miscellaneous changes for Matlab v5.

function [ngrps,grps,prob] = mstclust(crds,iter,critval,doplot)
  if (nargin < 2) iter = []; end;
  if (nargin < 3) critval = []; end;
  if (nargin < 4) doplot = []; end;

  if (isempty(iter))
    iter = 1000;
  end;
  if (isempty(critval))
    critval = 0.15;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  [N,P] = size(crds);
  if (P~=2)
    error('  Procedure is implemented only for 2 dimensions');
  end;
  nedge = N-1;

  [obs_edges,obs_edgelens,totlen] = mstree(crds); % Get observed MST
  [obs_edgelens,indx] = sort(obs_edgelens);       % Sort edge lengths
  obs_edgelens = flipud(obs_edgelens);
  obs_edges = flipud(obs_edges(indx,:));          % Carry along edges

  poly = hull(crds);                  % Get convex hull of data pts
  Nh = size(poly,1)-1;                % Number of points on hull
  Np = N - Nh;                        % Number of points in interior

  pts = zeros(N,2);                   % Allocate matrices
  prob = zeros(nedge,1);
  edgelen = zeros(iter,nedge);
%  treelen = zeros(iter,1);

  for it = 1:iter                     % Random trees
    pts = randcnvx(poly,N);             % Random pts within hull
    [e,randlen,t] = mstree(pts);        % MST of random pts
    randlen = randlen/t;                % Make relative to total tree length
    randlen = flipud(sort(randlen));    % Sort edge lengths
    edgelen(it,:) = randlen';           % Accum randomized edge lengths
%    treelen(it) = t;                    % Accum randomized tree lengths
  end;

%  medlen = median(treelen);           % Median randomized tree length
%  obs_edgelens = obs_edgelens*(medlen/totlen);  % Adjust observed tree length to median
  obs_edgelens = obs_edgelens / totlen;  % Make observed edges relative to total

  edgelen = sort(edgelen);
  prob = randprob(obs_edgelens,edgelen)';

  grps = ones(N,1);                   % Assume single group
  ngrps = 1;

  test = (prob<=(critval*ones(nedge,1))); % Check for significantly long edges
  if (sum(test)>0)
    [m,j] = min(test);
    if (j>1)
      ngrps = j;
      grps = mstgrp(crds,ngrps);
    end;
  end;

  if (doplot)
    figure(1); clf;
    hold on;
    plot(crds(:,1),crds(:,2),'ko');   % Plot grps with hulls
    putbnd(crds(:,1),crds(:,2));
    for i = 1:ngrps
      c = crds(grps==i,:);
      h = hull(c);
      plot(h(:,1),h(:,2),'k-');
    end;
    hold off;

    figure(2); clf;                   % Plot MST
    plot(crds(:,1),crds(:,2),'ko');
    putbnd(crds(:,1),crds(:,2));
    hold on;
    for i = 1:(N-1)
      t = obs_edges(i,1);
      u = obs_edges(i,2);
      x = [crds(t,1); crds(u,1)];
      y = [crds(t,2); crds(u,2)];
      plot(x,y,'k-');
      text(mean(x),mean(y),num2str(i));
    end;
    hold off;
  end;

  return;
