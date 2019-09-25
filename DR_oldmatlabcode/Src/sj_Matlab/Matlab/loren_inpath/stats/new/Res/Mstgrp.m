% MSTGRP: Find the k groups represented by the deletion of the k-1 longest
%         edges of a minimum spanning tree produced by MSTREE().
%
%     Syntax: grps = mstgrp(crds,ngrps,{doplot})
%
%         crds -   [n x p] matrix of point coordinates.
%         ngrps -  number of groups desired by deleting the ngrps-1 longest
%                    edges of the minimum spanning tree.
%         doplot - boolean flag, for the 2-dimensional case, indicating
%                     whether (=TRUE) or not (=FALSE) to plot the points and
%                     hulls [default = FALSE].
%         -----------------------------------------------------------------
%         grps -   [n x 1] vector of group identifiers; values are arbitrarily.
%

% RE Strauss, 12/31/96
%   9/7/99 - miscellaneous changes for Matlab v5.

function grps = mstgrp(crds,ngrps,doplot)
  if (nargin < 2) ngrps = []; end;
  if (nargin < 3) doplot = []; end;
    
  if (isempty(ngrps))
    error('MSTGRPS: Number of groups must be specified');
  end;
  if (isempty(doplot))
    doplot = 0;
  end;

  [edges,edgelen] = mstree(crds);

  [N,p] = size(crds);
  if (doplot & p~=2)
    doplot = FALSE;
  end;

  nedge = length(edgelen);
  n = nedge+1;
  grps = zeros(n,1);
  ndel = ngrps-1;
  v = ones(n,1);                        % Boolean list of vertices

%  % Heuristic: because the MST procedure tends to isolate outliers before
%  % identifying groups, do not delete a longest edge if it isolates a single 
%  % point; rather, skip it and proceed to the next longest edge.

  i = 0;
  e = edges(:);
  while (i < ndel)
    [m,j] = max(edgelen);
    edgelen(j) = 0;
    v1 = edges(j,1);
    v2 = edges(j,2);
%    if (sum(v1==e)>1 & sum(v2==e)>1)
%%      edges(j,1) = edges(j,2);
      edges(j,:) = [];
      edgelen(j) = [];
      e = edges(:);
      i = i+1;
%    end;
  end;

  edgelist = [edges; fliplr(edges)];    % Create sparse adjacency matrix
  adjacency = sparse(edgelist(:,1),edgelist(:,2),ones(size(edgelist,1),1));

  for k = 1:ngrps                       % Identify subset of vertices for each grp
    [m,j] = max(v);                       % First vertex of next grp
    visit = zeros(n,1);                   % Vector of vertices to be visited
    visit(j) = 1;
    visit = visitadj(j,visit,adjacency);  % Visit vertices
    i = find(visit);                      % Assign group identifiers
    leni = length(i);
    grps(i) = k * ones(leni,1);
    v(i) = zeros(leni,1);                 % Remove vertices from list
  end;

  if (doplot)                           % Plot 2D scatterplot of points and hulls
    clf;
%    plot(crds(:,1),crds(:,2),'ko');
    plotnum(crds(:,1),crds(:,2));
    hold on;
    for k = 1:ngrps
      X = crds(grps==k,:);
      h = hull(X);
      plot(h(:,1),h(:,2),'k');
    end;
    hold off;
  end;

  return;

