% GABCLUST: Constrained clustering on a Gabriel graph.
%           Produces plot of dendrogram if topology is complete; if incomplete,
%           a scatterplot is produced of either the two columns of the data matrix
%           (for p=2) or of the first two principal components of the correlation
%           matrix (for p>2).
%
%     Usage: [topology,nclust,grpid,support] = ...
%               gabclust(X,{labels},{iter},{suppress},{fontsize})
%
%         X =         [n x p] data matrix.
%         labels =    optional [n x q] matrix of group labels for dendrogram.
%         iter =      number of bootstrap iterations.
%         suppress =  optional flag; if present and 1 (=1), suppresses
%                       graphical dendrogram output.
%         fontsize =  optional font size for labels [default = 10].
%         -----------------------------------------------------------------------
%         topology =  [(n-1) x 4] matrix summarizing dendrogram topology:
%                       col 1 = 1st OTU/cluster being grouped at current step
%                       col 2 = 2nd OTU/cluster
%                       col 3 = ID of cluster being produced
%                       col 4 = distance at node
%         nclust =    estimated number of clusters (modes).
%                       If iter=0, nclust is the max number of clusters present at
%                         any step in the clustering process.
%                       If iter>0, nclust is ...
%         grpid =     [n x 1] vector of group-membership identifiers for all
%                       observations.
%         support =   [(n-2) x (n-1)] matrix, with one row for all but the base
%                       node, specifying group membership (support) at each node.
%

% RE Strauss, 1/17/99 (modified from Knnclust)
%   9/7/99 - changed plot colors for Matlab v5.

function [topology,nclust,grpid,support] = gabclust(X,labels,iter,suppress,fontsize)
  if (nargin < 2) labels = []; end;
  if (nargin < 3) iter = []; end;
  if (nargin < 4) suppress = []; end;
  if (nargin < 5) fontsize = []; end;

  [n,p] = size(X);                    % Size of data matrix

  get_nclust = 1;                     % Set output flags
  if (nargout < 2)
    get_nclust = 0;
  end;

  get_grpid = 1;
  if (nargout < 3)
    get_grpid = 0;
  end;

  if (nargout < 4)
    suprt = 0;
  else
    suprt = 1;
    support = zeros(n-2,n-1);
  end;

  if (~isempty(labels))               % Default input arguments
    if (size(labels,1)~=n)
      error('  Numbers of taxa and taxon labels do not match');
    end;
  end;

  if (isempty(iter))
    iter = 0;
  end;

  if (isempty(suppress))
    suppress = 0;
  end;

  topology = gab(X);                  % Gabriel clustering algorithm

  topo_compl = 1;                     % Complete-topology flag
  if (any(sum(topology')==0))
    topo_compl = 0;
  end;

  grpid = [];
  critval = [];

  if (get_nclust)                     % Estimate number of clusters
    [t,cl,nc] = topotips(topology);     % Evaluate topology

    if (~iter | ~topo_compl)            % Estimate by max number of clusters
%      nclust = max(nc);                   % Max number of clusters
nclust = nc(length(nc));
      if (get_grpid)
        i = max(find(nc==nclust));        % Find last link at which have max clusters
        grpid = cl(i,:);                  % Cluster membership

        i = find(grpid==0);               % Consecutively number terminal taxa
        if (~isempty(i))
          grpid(i) = [1:length(i)];
        end;

        u = uniquef(grpid);                % Convert cluster ids to 1,2,3...
        g = grpid;
        for iu = 1:length(u)
          i = find(g==u(iu));
          grpid(i) = iu*ones(length(i),1);
        end;
      end;

    else % if (iter)                    % Or by bootstrapping clustering nodes
knnd = [];
      critval = knnboot(X,3,knnd',iter,0.08);
critval

      i = min(find(topology(:,4)>=critval));  % Find actual nodes >= critical value
      if (isempty(i))                     % If none are significant,
        nclust = 1;                         % There's only one group
        grpid = ones(1,n);
      else                                % Else identify groups
        grpid = cl(i-1,:);

        i = find(grpid==0);                 % Consecutively number terminal taxa
        if (~isempty(i))
          grpid(i) = [1:length(i)];
        end;

        u = uniquef(grpid);
        nclust = length(u);

        if (get_grpid)
          g = grpid;                      % Convert cluster ids to 1,2,3...
            for iu = 1:length(u)
          i = find(g==u(iu));
            grpid(i) = iu*ones(length(i),1);
          end;
        end;
      end;
    end;
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

  if (~suppress)                      % Plot dendrogram
    if (topo_compl)
      figure;
      dendplot(topology,labels,fontsize);

      if (~isempty(critval))
        hold on;
        v = axis;
        plot([critval critval],v(3:4),'b');
        hold off;
      end;
    else
      if (isempty(grpid))
        [taxa,clust] = topotips(topology);
        [r,c] = size(clust);
        grpid = clust(r,:);
      end;

      if (p==2)
        figure;
        plotgrps(X(:,1),X(:,2),grpid');
        axis('equal');
      elseif (p>2)
        [loadings,percvar,scores] = pcacorr(X,2);
        figure;
        plotgrps(scores(:,1),scores(:,2),grpid');
        putxlab('PC1',percvar(1));
        putylab('PC2',percvar(2));
      else
        disp('GABCLUST: topology incomplete, no plot produced');
      end;
    end;
  end;

  return;
