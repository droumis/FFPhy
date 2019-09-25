% UPGMAF: Objective function for upgmaflx().
%
%     Usage: sse = upgmaf(beta,dist)
%
%         dist =  [n x n] symmetric distance matrix.
%         beta =  optional flexible-strategy parameter [default = 0].
%         -----------------------------------------------------------------------
%         sse =   squared sum of deviations of patristic distances from original
%                   distances.
%

% Belbin, L, DP Faith & GW Milligan. 1992. A comparison of two approaches to beta-
%   flexible clustering. Multivariate Behavioral Research 27:417-433.

% RE Strauss, 1/3/99

function sse = upgmaf(beta,dist)
  clstsize = ones(1,n);               % Number of elements in clusters/otus
  id = 1:n;                           % Cluster IDs
  topology = zeros(n-1,4);            % Output dendrogram-topology matrix

  plug = 10e6;
  dist = dist + eye(n)*plug;          % Replace diagonal with plugs

  for step = 1:(n-1)                  % Clustering steps
    min_dist = min(dist(:));            % Find minimum pairwise distance
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

    new_clstsize = clstsize(i) + clstsize(j);
    alpha_i = (1-beta) * clstsize(i) / new_clstsize;
    alpha_j = (1-beta) * clstsize(j) / new_clstsize;
    clstsize(i) = new_clstsize;

    for k = 1:n                         % For all other clusters/OTUs,
      if (k~=i & k~=j)                %   adjust distances to new cluster
        dist(k,i) = alpha_i*dist(k,i) + alpha_j*dist(k,j) + beta*min_dist;
        dist(i,k) = alpha_i*dist(i,k) + alpha_j*dist(j,k) + beta*min_dist;
        dist(k,j) = plug;
        dist(j,k) = plug;
      end;
    end;
  end; % for step


  return;
