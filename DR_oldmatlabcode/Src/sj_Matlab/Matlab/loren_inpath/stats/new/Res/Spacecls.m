% SPACECLS: K-means clustering for 2D spatial distribution.  Varies the numbers
%           of clusters from 2 to a specified maximum and evaluates the
%           clustering solution with respect with randomly distributed data
%           within the observed convex hull of the data.  Minimizes the within-
%           group sum of squared distances from group centroids.
%
%     Usage: [sse,prob] = spacecls(crds,maxk,restarts)
%
%           crds = [n x 2] matrix of spatial coordinates for n individuals
%           maxk = maximum number of clusters evaluated
%           restarts = number of randomized restarts per cluster analysis
%                       [default=0]
%
%           sse =  [maxk x 1] vector of SSE's for 1:maxk clusters
%           prob = [maxk x 1] vector of corresponding probabilities with
%                   respect to k-means clustering of randomized spatial
%                   distributions
%

% RE Strauss, 4/27/97

function [sse,prob] = spacecls(crds,maxk,restarts)
  tic;
  if (nargin < 3)
    restarts = 0;
  end;

  [n,p] = size(crds);
  randcrds = crds;                    % Matrix for random coordinates
  ndistrib = 50;                      % Number of randomizations for probability

  sse = zeros(maxk,1);
  prob = zeros(maxk,1);

  h = hull(crds);                     % Convex hull
  nhull = length(h);                  % Number of hull points
  mincrds = min(h);                   % Extents of point distribution
  maxcrds = max(h);
  rangecrds = maxcrds - mincrds;

  centr = mean(crds);                 % Centroid
  dist = eucl(crds,centr);            % Euclidean distances from centroid
  sse(1) = dist' * dist;              % Total sum-of-squares

  hullarea = polyarea(h);             % Area of hull
  radius = sqrt(hullarea/pi);         % Radius of circle of same area

%  radius = mean(eucl(h,centr));

k=1
sse

plot(crds(:,1),crds(:,2),'ok');
hold on;
plot(h(:,1),h(:,2),'k');
plot(centr(1),centr(2),'+k');
axis('equal');
axis('square');
hold off;
%pause;

  for r = 1:ndistrib                  % Probability of total sse
    randcrds = randcirc(n,radius,centr);


%plot(randcrds(:,1),randcrds(:,2),'ok');
%hold on;
%plot(h(:,1),h(:,2),'k');
%plot(centr(1),centr(2),'+k');
%axis('equal');
%axis('square');
%hold off;
%pause;

    centr = mean(randcrds);
    dist = eucl(randcrds,centr);
    s = dist' * dist;
%s
    if (s < sse(1))
      prob(1) = prob(1) + 1/ndistrib;
    end;
  end;
sse
prob

  for k = 2:maxk
k
    [centr,clst,s] = kmeans(crds,k,restarts);
    sse(k) = s;
s

    for r = 1:ndistrib                  % Probability of total sse
      randcrds = randcirc(n,radius,mean(crds));
      [centr,clst,s] = kmeans(randcrds,k);

plot(randcrds(:,1),randcrds(:,2),'ok');
hold on;
plot(h(:,1),h(:,2),'r');
uc = uniquef(clst);
for uci = 1:length(uc)
  ucv = uc(uci);
  ucindx = find(clst==ucv);
  ucx = randcrds(ucindx,:);
  uch = hull(ucx);
  plot(uch(:,1),uch(:,2),'k');
end;
hold off;
%pause;

s
      if (s < sse(k))
        prob(k) = prob(k) + 1/ndistrib;
      end;
    end;
sse
prob
  end;

  minutes = toc/60
  return;
