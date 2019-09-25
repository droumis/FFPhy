% KNNPLOT: Produces either a scatterplot or a dendrogram for a topology produced
%           by kNN cluster analysis, depending on whether the topology is complete.
%
%     Usage: knnplot(crds,topology,labels,fontsize)
%
%         crds =      [n x p] data matrix.
%         topology =  [(n-1) x 4] matrix summarizing dendrogram topology:
%                       col 1 = 1st OTU/cluster being grouped at current step
%                       col 2 = 2nd OTU/cluster
%                       col 3 = ID of cluster being produced
%                       col 4 = distance at node
%         labels =    optional [n x q] matrix of group labels for dendrogram.
%         fontsize =  optional font size for labels [default = 10].
%

% RE Strauss, 1/26/99, extracted from knnclust.m

function knnplot(crds,topology,labels,fontsize)
  if (nargin < 3)
    labels = [];
  end;
  if (nargin < 4)
    fontsize = [];
  end;

  [n,p] = size(crds);

  topo_compl = 1;                     % Complete-topology flag
  if (any(sum(topology')==0))
    topo_compl = 0;
  end;

  if (topo_compl)
    figure;
    dendplot(topology,labels,fontsize);

  else
    [taxa,clust] = topotips(topology);
    [r,c] = size(clust);
    grpid = clust(r,:);

    if (p==2)
      figure;
      plotgrps(crds(:,1),crds(:,2),grpid');
      axis('equal');
    else
      [loadings,percvar,scores] = pcacorr(crds,2);
      figure;
      plotgrps(scores(:,1),scores(:,2),grpid');
      putxlab('PC1',percvar(1));
      putylab('PC2',percvar(2));
    end;
  end;

  return;

