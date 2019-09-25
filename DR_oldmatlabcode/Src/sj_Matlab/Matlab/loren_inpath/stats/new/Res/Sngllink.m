% SNGLLINK: Single-linkage hierarchical cluster analysis of a distance
%         matrix.  Produces plot of dendrogram.
%
%     Usage: [topology,support] = sngllink(dist,{labels},{suppress},{fontsize})
%
%         dist =     [n x n] symmetric distance matrix.
%         labels =   optional [n x q] matrix of group labels for dendrogram.
%         suppress = optional flag; if present and true (=1), suppresses
%                       graphical dendrogram output.
%         fontsize = optional font size for labels [default = 10].
%         -----------------------------------------------------------------------
%         topology = [(n-1) x 4] matrix summarizing dendrogram topology:
%                       col 1 = 1st OTU/cluster being grouped at current step
%                       col 2 = 2nd OTU/cluster
%                       col 3 = ID of cluster being produced
%                       col 4 = distance at node
%         support =  [(n-2) x (n-1)] matrix, with one row for all but the base 
%                       node, specifying group membership (support) at each node.
%

% RE Strauss, modified from upgma() 1/2/99
%   9/7/99 - miscellaneous changes for Matlab v5.

function [topology,support] = sngllink(dist,labels,suppress,fontsize)
  if (nargin < 2) labels = []; end;
  if (nargin < 3) suppress = []; end;
  if (nargin < 4) fontsize = []; end;

  [n,p] = size(dist);
  if (n~=p | any(diag(dist)))
    dist
    error('  Input matrix is not a distance matrix');
  end;

  if (~isempty(labels))
    if (size(labels,1)~=n)
      error('  Numbers of taxa and taxon labels do not match');
    end;
  end;

  if (isempty(suppress))
    suppress = 0;
  end;

  if (nargout<2)
    suprt = 0;
  else
    suprt = 1;
    support = zeros(n-2,n-1);
  end;

  id = 1:n;                           % Cluster IDs
  topology = zeros(n-1,4);            % Output dendrogram-topology matrix

  plug = max(max(dist)) + 1e6;
  dist = dist + eye(n)*plug;          % Replace diagonal with plugs

  for step = 1:(n-1)                  % Clustering steps
    min_dist = min(dist(:));            % Find minimum pairwise distance
if (~finite(min_dist))
  step
  min_dist
end;
    [ii,jj] = find(dist==min_dist);     % Find location of minimum
    k = 1;                              % Use first identified minimum
    while (ii(k)>jj(k))                 %   for which i<j
      k = k+1;
    end;
    i = ii(k);
    j = jj(k);
    if (id(i)<id(j))
      topology(step,:) = [id(i) id(j) n+step min_dist];
    else
      topology(step,:) = [id(j) id(i) n+step min_dist];
    end;
    id(i) = n+step;
    dist(i,j) = plug;
    dist(j,i) = plug;

    for h = 1:n                         % For all other clusters/OTUs,
      if (h~=i & h~=j)                  %   adjust distances to new cluster
        dhi = dist(h,i);
        dhj = dist(h,j);

        dist(h,i) = min([dhi dhj]);
        dist(i,h) = dist(h,i);
        dist(h,j) = Inf;
        dist(j,h) = Inf;
      end;
    end;
  end; % for step

  if (~suppress)                      % Plot dendrogram
    dendplot(topology,labels,fontsize);
  end;

  if (suprt)                          % Specify group membership at nodes
    for node = 1:size(topology,1)-1
      pos = 0;
      t1 = topology(node,1);
      t2 = topology(node,2);

      pos = pos+1;
      if (t1 <= n)
        support(node,pos) = t1;
      else
        mem = find(support(t1-n,:)>0);
        support(node,pos:(pos+length(mem)-1)) = support(t1-n,mem);
        pos = pos+length(mem)-1;
      end;
    
      pos = pos+1;
      if (t2 <= n)
        support(node,pos) = t2;
      else
        mem = find(support(t2-n,:)>0);
        support(node,pos:(pos+length(mem)-1)) = support(t2-n,mem);
        pos = pos+length(mem)-1;
      end;

      len = length(find(support(node,:)>0));
      support(node,1:len) = sort(support(node,1:len));
    end;
  end;

  return;
