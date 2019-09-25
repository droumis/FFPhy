% KNNBOOT:  Bootstrap the kth nearest-neighbor procedure to estimate the number
%           of clusters and corresponding group-membership of terminal taxa.
%           Use the bootstrap method of Wong (1985), but modified to use
%           dendrogram-node distance as a test statistic.
%           Called by knnclust().
%
%     Usage: critval = knnboot(X,k,knnd,iter,alpha)
%

function critval = knnboot(X,k,knnd,iter,alpha)
  [n,p] = size(X);

  if (alpha > 1)                        % Convert alpha to 1-(tail proportion)
    alpha = 0.01*alpha;
  end;
  if (alpha < 0.5)
    alpha = 1-alpha;
  end;

  maxdistrib = (n-2)*iter;
  distrib = zeros(maxdistrib,1);
  ib = 1;

  for it = 1:iter                       % For each bootstrap iteration,
    r = ceil(n*rand(n,1));                % Random indices
    while (std(r)==0)                     % If all indices identical,
      r = ceil(n*rand(n,1));              %   resample
    end;
    knndb = knnd(r)*ones(1,p);              % Rearrange into bootstrapped sample
    u = rand(length(knnd),p)*2-1;           % Uniform-random on interval [-1,1]
    Xb = X(r,:) + knndb.*u;                 % As per Wong (1985)

    topo = knn(Xb,k);                     % KNN topology
    d = topo((topo(:,4)>0),4);            % Link distances

    ie = ib+length(d)-1;                  % Indices into distrib
    distrib(ib:ie) = d;                   % Stash branch lengths into distribution
    ib = ie+1;
  end;

  if (maxdistrib > ib)
    distrib(ib:maxdistrib) = [];        % Delete unused distribution buffer
  end;
  distrib = sort(distrib);              % Sort accumulated br-len distribution
  lend = length(distrib);

  c = lend*alpha;                       % Critical index + fraction
  cu = floor(c);
  cl = ceil(c);
  fr = 0;
  if (cu > cl)
    fr = (c-cu)/(cl-cu);
  end;

  du = distrib(cu);                     % Linearly interpolate critical value
  dl = distrib(cl);
  critval = du + fr*(dl-du);

  return;
