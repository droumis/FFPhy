% KNN: Kth nearest-neighbor clustering procedure, using method of ...  This function
%       produces the topology matrix from the input data matrix and value of k,
%       and does no error-checking.
%       See knnclust() for error-checking, plots, bootstrapping, and supporting
%       information.
%
%     Usage: [topology,knnd] = knn(X,k)
%
%         X =         [n x p] data matrix.
%         k =         smoothing parameter, the number of nearest neighbors included
%                       in the distance function.
%         -----------------------------------------------------------------------------
%         topology =  [(n-1) x 4] matrix summarizing dendrogram topology:
%                       col 1 = 1st OTU/cluster being grouped at current step
%                       col 2 = 2nd OTU/cluster
%                       col 3 = ID of cluster being produced
%                       col 4 = distance at node
%         knnd =      [n x 1] vector of kth nearest-neighbor distances of observations.
%

% RE Strauss, 1/7/99

function [topology,knnd] = knn(X,k)
  [n,p] = size(X);

  de = eucl(X);                          % Euclidean distances among points
  de = putdiag(de,Inf);

  knnd = sort(de);                       % Sort distances to points
  knnd = knnd(k,:);                      % Extract kth nearest-neighbor dists

  dist = zeros(n,n);                     % Initialize distance matrix
  dist = putdiag(dist,Inf);

  for i = 1:(n-1)
    for j = (i+1):n
      distij = de(i,j);
      if (distij<=knnd(i) | distij<=knnd(j))    % If pts i & j are neighbors,
        dist(i,j) = 0.5*(knnd(i)+knnd(j));      %   dist is mean of their NNDs
        dist(j,i) = dist(i,j);
      else                                      % Else dist is infinity
        dist(i,j) = Inf;
        dist(j,i) = Inf;
      end;
    end;
  end;

  id = 1:n;                           % Cluster IDs
  topology = zeros(n-1,4);            % Output dendrogram-topology matrix

  for step = 1:(n-1)                  % Single-linkage clustering steps
    min_dist = min(dist(:));            % Find minimum pairwise distance

    if (finite(min_dist))               % If a min distance exists,
      [ii,jj] = find(dist==min_dist);     % Find location of minimum
      h = 1;                              % Use first identified minimum
      while (ii(h)>jj(h))                 %   for which i<j
        h = h+1;
      end;
      i = ii(h);
      j = jj(h);
      dij = min_dist;

      if (id(i)<id(j))
        topology(step,:) = [id(i) id(j) n+step dij];
      else
        topology(step,:) = [id(j) id(i) n+step dij];
      end;
      id(i) = n+step;
      dist(i,j) = Inf;
      dist(j,i) = Inf;

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
    end;
  end; % for step

  return;
