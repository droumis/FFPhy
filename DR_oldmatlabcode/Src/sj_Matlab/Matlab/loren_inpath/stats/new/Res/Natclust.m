% NATCLUST: Finds "natural" clusters using the algorithm of Carmichael et al (1968).
%
%     Usage: [ngrps,grps] = natclust(crds,{ngrps},{doplot})
%
%           crds -    [n x p] matrix of point coordinates.
%           ngrps -   optional number of groups desired [default = all].
%           doplot -  optional boolean value indicating whether (1) or not 
%                       (0, default) to plot points and hulls.
%           ------------------------------------------------------------------------
%           ngrps -   [1 x nres] vector of numbers of groups recognized at each 
%                       resolution level.
%           grps -    [n x nres] matrix of assigned group memberships for all points
%                       at each resolution level.  Grp-id = 0 for unassigned 
%                       points.
%

% Carmichael JW, JA George & RS Julius. 1968. Finding natural clusters.
%   Syst. Zool. 17:144-150.

% RE Strauss, 11/12/98

function [ngrps,grps] = natclust(crds,ngrps,doplot)
  if (nargin < 2) ngrps = []; end;
  if (nargin < 3) doplot = []; end;

  target_ngrps = [];
  if (~isempty(ngrps))
    target_ngrps = ngrps;
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  R = 15;                                 % Number of resolution levels
  maxres = 0.20;                          % Percentiles of max & min resolution levels
  minres = 0.90;

  [N,P] = size(crds);
  grps =  zeros(N,R);                     % Allocate output matrix

  dist = eucl(crds);                      % Euclidean distances among points
  [dist,pt1,pt2] = trilow(dist);          % Lower triangular distance matrix
  [dist,pt1,pt2] = sortmat(dist,pt1,pt2); % Sort dists low to high, carry along pt ids
  ndist = length(dist);

  hyperdiag = eucl(min(crds),max(crds));  % Length of hyperdiagonal of axis ranges
  dist = dist./hyperdiag;                 % Normalized distances

  maxres = max([floor(maxres*ndist),2]);    % Set resolution levels
  minres = min([ceil(minres*ndist),ndist-1]);
  res1 = linspace(dist(minres),dist(maxres),R)';
  max_res1 = res1(1);
  res1_range =  max_res1-res1(R);
  res2_lowval = max_res1 - 0.5*res1_range;
  res2 = linspace(max_res1,res2_lowval,R)';

  for res_level = 1:R                     % Cycle thru resolution levels
    g = zeros(N,1);                         % Initialize group-membership vector
    
    d = dist;                               % Initialize distance & pt-id vectors
    p1 = pt1;
    p2 = pt2;
    r1 = res1(res_level);                 % Current resolution values
    r2 = res2(res_level);

    grpid = 0;
    while (d(1) < r1)                       % Assign points to clusters
      grpid = grpid+1;                        % Next cluster
      curgrp = [p1(1) p2(1)];                 % Grab 2 closest points
      avglink = d(1);                         % Average link
      slink = d(1);                           % Single link
      d(1) = [];                              % Delete pts from distance & pt-id vectors
      p1(1) = [];
      p2(1) = [];

      newpoint = 0;
      if (~isempty(d))
        newpoint = 1;                         % Current cluster
      end;
      while (newpoint)                        % Find pt closest to cluster
        i = min(find(isin(p1,curgrp)));
        j = min(find(isin(p2,curgrp)));
        if (isempty(i))
          cpi = p1(j);
          if (isempty(j))
            newpoint = 0;
          end;
        elseif (isempty(j))
          cpi = p2(i);
        elseif (i<j)                          % Crd index of closest point
          cpi = p2(i);
        else
          cpi = p1(j);
        end;
        if (newpoint)
          dist_to_new = eucl(crds(curgrp,:),crds(cpi,:))./hyperdiag; % Dists from new pt to old pts

          new_avglink = mean(dist_to_new);      % Average-link criterion
          drop = avglink - new_avglink;
          diff = new_avglink - drop;
          avglink = new_avglink;
          if (diff > r1)
            newpoint = 0;
          end;

          new_slink = min(dist_to_new);         % Single-link criterion
          drop = mean(slink) - new_slink;       
          diff = new_slink - drop;
          slink = [slink new_slink];
          if (diff > r2)
            newpoint = 0;
          end;

          new_complink = max(dist_to_new);      % Ratio of complete link to single link
          linkratio = new_complink / new_slink;
          if (linkratio < r2)
            newpoint = 0;
          end;
        end;

        if (newpoint)
          curgrp = [curgrp cpi];              % Add point to cluster
          i = find(isin(p1,curgrp) & isin(p2,curgrp));
          if (~isempty(i))
            d(i) = [];                        % Delete pt from distance & pt-id vectors
            p1(i) = [];
            p2(i) = [];
          end;
        end;
      end;

      g(curgrp) = grpid * ones(length(curgrp),1); % Stash group labels
      i = find(isin(p1,curgrp) | isin(p2,curgrp));
      if (~isempty(i))
        d(i) = [];                            % Delete pts from distance & pt-id vectors
        p1(i) = [];
        p2(i) = [];
      end;

      if (isempty(d))
        d = 1e6;
      end;
    end;

    i = find(g==0);                         % Locate unassigned points
    nextgrp = max(g)+1;                     % Next group identifier
    g(i) = nextgrp:(nextgrp+length(i)-1);   % Put unassigned points into their own grps

    grps(:,res_level) = g;                  % Stash group-membership vector
  end;

  g = grps;                                 % Find set of unique solutions
  n = ngrps;

  j = 1;
  u = uniquef(g(:,1));
  while (length(u)==1)
    j = j+1;
    u = uniquef(g(:,j));
  end;
  grps = g(:,j);

  lastcomp = 1;
  for i = j+1:R
    if (abs(sum(g(:,i)-grps(:,lastcomp)))>0)
      grps = [grps g(:,i)];
      lastcomp = lastcomp+1;
    end;
  end;

  [n,r] = size(grps);

  if (doplot)                               % Optional plots of groups
    close all;
    for i = 1:r
      figure;
      plotgrps(crds(:,1),crds(:,2),grps(:,i),[],[1 1 0]);
    end;
  end;

  ngrps = zeros(1,r);                       % Find numbers of clusters
  for i = 1:r                               % Convert singlet labels to zero
    [u,f] = uniquef(grps(:,i));
    if (any(f==1))
      for j = 1:length(u)
        if (f(j)==1)
          k = find(grps(:,i)==u(j));
          grps(k,i) = zeros(length(k),1);
        end;
      end;
      ngrps(i) = length(uniquef(grps(:,i)))-1;
    else
      ngrps(i) = length(u);
    end;    
  end;

  if (~isempty(target_ngrps))               % Restrict to specified number of groups
    i = find(ngrps==target_ngrps);
    ngrps = ngrps(i);
    grps = grps(:,i);
  end;
  
  return;

