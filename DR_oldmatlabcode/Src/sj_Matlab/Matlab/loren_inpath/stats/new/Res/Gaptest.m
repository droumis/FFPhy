% GAPTEST:  A pairwise test of the significance of the gaps among the convex hulls
%           of two or more clusters of points (Rasson & Kubushishi 1994), under
%           the null hypothesis that two clusters are random samples from a 
%           single cluster.
%
%     Usage: [prob,hullarea,gapspace,npts,lambda] = gaptest(g,X)
%
%           g = grouping vector (length n) for k clusters of points.
%           X = [n x 2] matrix of two-dimensional point coordinates.
%           -----------------------------------------------------------------------
%           prob =      [k x k] matrix of pairwise levels of statistical 
%                         significance.
%           hullarea =  corresponding matrix of convex-hull areas: within-cluster
%                         estimates on the diagonal, pairwise pooled-cluster 
%                         on the off-diagonals.
%           gapspace =  corresponding matrix of gap-space estimates, m'.
%           npts =      corresponding matrix of sample sizes.
%           lambda =    corresponding estimates of the stationary Poisson point-
%                         process intensity parameter.
%

% Rasson, J.-P. & T. Kubushishi. 1994. The gap test: an optimal method for
%   determining the number of natural clusters in cluster analysis.
%   In: E. Diday, Y. Lechevallier, M. Schader, P. Bertrand & B. Burtschy (eds.),
%   New Approaches in Classification and Data Analysis, pp. 186-193.
%   Springer-Verlag.

% RE Strauss, 1/3/99

function [prob,hullarea,gapspace,npts,lambda] = gaptest(g,X)
  [r,c] = size(g);
  if (min([r,c])>1)
    error('GAPTEST: group-identification vector must be a vector');
  end;
  ng = max([r,c]);

  clustid = uniquef(g);
  nclust = length(clustid);

  [n,p] = size(X);
  if (p~=2)
    error('GAPTEST: two-dimensional points only');
  end;
  if (n~=ng)  
    error('GAPTEST: input matrices not compatible');
  end;

  hullarea = zeros(nclust,nclust);        % Allocate output matrices
  gapspace = zeros(nclust,nclust);

  hullpts = [];                           % Matrices of hull points and identifiers
  hid = [];
  npts = zeros(nclust,nclust);            % Numbers of points per cluster                           
    
  for ik = 1:nclust                       % Estimate within-cluster parameters
    i = find(g==clustid(ik));               % Isolate points for cluster
    npts(ik,ik) = length(i);
    h = hull(X(i,:));
    hullpts = [hullpts; h];                 % Append hull points to list
    hid = [hid; ik*ones(size(h,1),1)];
    hullarea(ik,ik) = polyarea(h);          % Area of convex hull
  end;

  for i = 1:(nclust-1)                    % Estimate pooled-cluster parameters
    h1 = hullpts(hid==i,:);

    for j = (i+1):nclust
      h2 = hullpts(hid==j,:);

      npts(i,j) = npts(i,i)+npts(j,j);
      npts(j,i) = npts(i,j);

      h = hull([h1;h2]);
      hullarea(i,j) = polyarea(h);          % Area of pooled hulls
      hullarea(j,i) = hullarea(i,j);
                                            % Gap-space statistic m'
      gapspace(i,j) = hullarea(i,j) - hullarea(i,i) - hullarea(j,j);
      gapspace(j,i) = gapspace(i,j);
    end;
  end;

  lambda = npts./hullarea;                % Poisson lambda estimate
  prob = exp(-lambda.*gapspace);

  return;
