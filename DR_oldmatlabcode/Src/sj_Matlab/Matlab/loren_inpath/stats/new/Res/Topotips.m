% TOPOTIPS: Finds the terminal taxa within clusters from a topology matrix.  
%           The topology matrix may be incomplete, in the sense of omitting 
%           some higher-order clusters.
%
%     Usage: [taxa,clusters,nclust,ntaxa,clust_id] = topotips(topology)
%
%         topology =  [(n-1) x 4] matrix summarizing dendrogram topology:
%                       col 1 = 1st OTU/cluster being grouped at current step
%                       col 2 = 2nd OTU/cluster
%                       col 3 = ID of cluster being produced
%                       col 4 = distance at node
%         --------------------------------------------------------------------
%         taxa =      [(n-1) x n] matrix of taxa (rows) in each cluster 
%                       (or fewer rows for an incomplete topology matrix).
%         clusters =  matrix of cluster identifier for each taxon at each clustering 
%                       step.
%         nclust =    number of clusters (excluding terminal taxa) at each 
%                       clustering step.
%         ntaxa =     number of terminal taxa per cluster.
%         clust_id =  column vector of cluster identifiers.
%

% RE Strauss, 1/2/99

function [taxa,clusters,nclust,ntaxa,clust_id] = topotips(topology)
  [r,p] = size(topology);
  ntaxa = topology(1,3)-1;                % Number of terminal taxa
  nclust = max(topology(:,3))-ntaxa;      % Number of clusters

  clust_id = [(ntaxa+1):(nclust+ntaxa)]'; % Allocate output matrices
  taxa =     zeros(nclust,ntaxa);
  clusters = zeros(nclust,ntaxa);
  nclust =   zeros(nclust,1);

  for ti = 1:r                            % Cycle thru topology, identifying taxa
    nt = 1;
    for g = 1:2                             % Two groups comprising cluster
      a = topology(ti,g);
      if (a <= ntaxa)                         % Group is terminal taxon
        taxa(ti,nt) = a;
        nt = nt+1;
      else                                    % Else group is previous cluster
        a = a-ntaxa;
        i = find(taxa(a,:)>0);
        taxa(ti,nt:(nt+length(i)-1)) = taxa(a,i);
        nt = nt+length(i);
      end;
    end;
  end;

  for ti = 1:r                            % Cycle thru topology, identifying clusters
    if (ti>1)
      clusters(ti,:) = clusters(ti-1,:);
    end;
    t = find(taxa(ti,:)>0);
    for i = 1:length(t)
      tt = taxa(ti,i);
      clusters(ti,tt) = clust_id(ti);
    end;
    c = find(clusters(ti,:)>0);
    nclust(ti) = length(uniquef(clusters(ti,c)));
  end;

  for ti = 1:r                            % Cycle thru again, inserting single taxa
    c = find(clusters(ti,:)==0);
    clusters(ti,c) = c;
  end;

  coltot = sum(taxa>0);                   % Compress blank cols & rows
  i = find(coltot==0);
  taxa(:,i) = [];

  rowtot = sum((taxa>0)')';
  i = find(rowtot==0);
  taxa(i,:) = [];

  ntaxa = sum((taxa>0)')';                % Number of taxa per cluster
  return;

