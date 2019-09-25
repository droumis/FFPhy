% MSTREE: Minimum spanning tree
%
%     Syntax: [edges,edgelen,totlen] = mstree(crds,{doplot})
%
%         crds -    [n x p] matrix of point coordinates.
%         doplot -  boolean flag indicating whether (=1) or not (=0)
%                     to plot the points and minimum spanning tree
%                     [default = 0].
%         -----------------------------------------------------------------
%         edges -   [(n-1) x 2] list of points defining edges, in arbitrary
%                     sequence.
%         edgelen - [(n-1) x 1] vector of edge lengths corresponding to
%                     'edges'.
%         totlen -  total length of tree.
%

% RE Strauss, 7/9/96
%   2/9/99 - produces plot only for 2D coordinates, but will find the tree
%              for any number of dimensions.
%   9/7/99 - changed plot colors for Matlab v5.

function [edges,edgelen,totlen] = mstree(crds,doplot)
  if (nargin < 2) doplot = []; end;

  [n,p] = size(crds);

  if (isempty(doplot) | p~=2)
    doplot = 0;
  end;

  totlen = 0;
  edges = zeros(n-1,2);
  edgelen = zeros(n-1,1);

  dist = eucl(crds);            % Pairwise euclidean distances
  highval = max(max(dist))+1;


  e1 = 1*ones(n-1,1);           % Initialize potential-edge list
  e2 = [2:n]';
  ed =  dist(2:n,1);

  for edge = 1:(n-1)            % Find shortest edges
    [mindist,i] = min(ed);        % New edge
    t = e1(i);
    u = e2(i);
    totlen = totlen + mindist;
    if (t<u)                      % Add to actual-edge list
      edges(edge,:) = [t u];
    else
      edges(edge,:) = [u t];
    end;
    edgelen(edge) = mindist;

    if (edge < n-1)
      i = find(e2==u);              % Remove new vertex from
      e1(i) = 0;                    %   potential-edge list
      e2(i) = 0;
      ed(i) = highval;

      indx = find(e1>0);
      for i = 1:length(indx)        % Update potential-edge list
        j = indx(i);
        t = e1(j);
        v = e2(j);
        if (dist(u,v) < dist(t,v))
          e1(j) = u;
          ed(j) = dist(u,v);
        end;
      end;
    end;
  end;

  if (doplot)
    figure;
    plot(crds(:,1),crds(:,2),'ok');
    putbnd(crds(:,1),crds(:,2));
    hold on;
    for i = 1:(n-1)
      t = edges(i,1);
      u = edges(i,2);
      x = [crds(t,1); crds(u,1)];
      y = [crds(t,2); crds(u,2)];
      plot(x,y,'-k');
    end;
    hold off;
  end;

  return;
