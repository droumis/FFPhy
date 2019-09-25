% GAB: Constrained Gabriel clustering procedure.  This function
%       produces the topology matrix from the input data matrix,
%       and does no error-checking.
%       See gabclust() for error-checking, plots, bootstrapping, and supporting
%       information.
%
%     Usage: topology = gab(X)
%
%         X =         [n x p] data matrix.
%         -----------------------------------------------------------------------------
%         topology =  [(n-1) x 4] matrix summarizing dendrogram topology:
%                       col 1 = 1st OTU/cluster being grouped at current step
%                       col 2 = 2nd OTU/cluster
%                       col 3 = ID of cluster being produced
%                       col 4 = distance at node
%

% RE Strauss, 1/7/99

function topology = gab(X)
  [n,p] = size(X);

  [connect,dist] = gabriel(X,1);      % Get distances on Gabriel graph

  [i,j] = find(~connect);             % Set non-connected distances to Inf
  for k=1:length(i)
    dist(i(k),j(k)) = Inf;
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
