% KNNCLUST: Finds the best partition of point coordinates by the Kth nearest-neighbor
%           clustering procedure (Strauss 2001: 
%           <http://www.biol.ttu.edu/Strauss/Pubs/Papers/01AnimBehav.pdf>).
%
%     Usage: [nclust,grpid,adr,t,tprob,U,Uprob,k] = knnclust(crds,{doplots},{k})
%
%         crds =    [n x p] data matrix.
%         doplots = optional flag; if present and true (=1), produces graphical 
%                     dendrogram output [default = 0 = false].  For 2D data (p=2),
%                     scatterplots of the original data are produced; for >2D data,
%                     scatterplots are of principal components of the data. 
%         k =       optional smoothing parameter, the number of nearest neighbors 
%                     included in the distance function (range 1-n).  If not 
%                     provided, the value is varied from 1 to 5 to estimate the 
%                     number of groups present.
%         -------------------------------------------------------------------------
%         nclust =  estimated number of clusters (modes).
%         grpid =   [n x 1] vector of group-membership identifiers for all
%                     observations.
%         adr =     ADR statistic value for best cluster partition.
%         t =       corresponding t-statistic value for ADR.
%         tprob =   corresponding probability for t-statistic.
%         U =       corresponding U-statistic value for ADR.
%         Uprob =   corresponding probability for U-statistic.
%         k =       input value of k (if provided) or value corresponding to
%                     nclust (if not provided).
%

% RE Strauss, 12/24/98
%   10/22/02 - updated documentation and graphical output;
%              replace crds by PC scores for >2D.

% Note: Uses the wrong null distribution; must be based on predictions of this function from
%       random data.

function [nclust,grpid,adr,t,tprob,U,Uprob,k] = knnclust(crds,doplots,k)
  if (nargin < 2) doplots = []; end;
  if (nargin < 3) k = []; end;
  
  doprob = 0;
  if (nargout>3)
    doprob = 1;
  end;

  if (isempty(doplots))
    doplots = 0;
  end;

  [n,p] = size(crds);                   

  if (isempty(k))                                 % k not provided
    kmin = 1;
    kmax = min([5 n-1]);
  else                                            % k provided
    kmin = k;
    kmax = k;
  end;
  
  z = zeros(1,kmax-kmin+1);                       % Initialize accumulation matrices
  k_nclust = z;           
  k_adr = z;
  k_t = z;
  k_tprob = z;
  k_U = z;
  k_Uprob = z;
  k_k = z;
  k_grpid = zeros(n,kmax-kmin+1);

  compl_k = 0;
  topo_compl = 0;
  ik = 0;

  for k = kmin:kmax
    ik = ik+1;
    [topology,knnd] = knn(crds,k);                % KNN clustering algorithm
    if (all(rowsum(topology)>0))
      topo_compl = 1;
    end;

    if (topo_compl & compl_k==0)
      compl_k = k;
    end;

    [t,cl,nc] = topotips(topology);               % Evaluate topology
    if (topo_compl)
      cl(length(nc),:) = [];
    end;
    cl = cl';
    
    if (doprob)                                   % Estimate number of clusters
      [adr,nc,tprob,Uprob,df,t,U] = clustnum(crds,cl);    
      [minprob,i] = min(Uprob);
      adr = adr(i);
      nclust = nc(i);
      grpid = cl(:,i);
      t = t(i);
      tprob = tprob(i);
      U = U(i);
      Uprob = Uprob(i);
    else
      [adr,nc] = clustnum(crds,cl);
      [maxadr,i] = max(adr);
      adr = adr(i);
      nclust = nc(i);
      grpid = cl(:,i);
      t = NaN;
      tprob = NaN;
      U = NaN;
      Uprob = NaN;
    end;

    u = uniquef(grpid);                           % Convert cluster ids to 1,2,3...
    grpid = replace(grpid,u,1:length(u));

    k_nclust(ik) = nclust;
    k_adr(ik) = adr;
    k_t(ik) = t;
    k_tprob(ik) = tprob;
    k_U(ik) = U;
    k_Uprob(ik) = Uprob;
    k_k(ik) = k;
    k_grpid(:,ik) = grpid;
  end;

  [minprob,i] = min(k_Uprob);
  adr = k_adr(i);
  t = k_t(i);
  tprob = k_tprob(i);
  U = k_U(i);
  Uprob = k_Uprob(i);
  nclust = k_nclust(i);
  grpid = k_grpid(:,i);

  if (doplots)                                    % Graphical output
    if (p>2)                                        % For >2D, replace crds by PC scores
      [loadings,percvar,crds] = pcacorr(crds,2);
    end;
      
    figure;
    plotlabl(crds,tostr(1:n));                      % Scatterplot of data
    putbnds(crds);
    axis('equal');
    if (p>2)
      putxlab('PC1');
      putylab('PC2');
    end;
    
    figure;
    plotgrps(crds(:,1),crds(:,2),grpid,tostr(uniquef(grpid)));
    axis('equal');
    puttick('off','off');
    if (p>2)
      putxlab('PC1');
      putylab('PC2');
    end;
  end;

  return;
