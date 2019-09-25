% MahalDist: Finds the Mahalanobis distances of one or more points from a group
%            centroid, based on the estimated covariance matrix of the group.
%
%     D2 = mahaldist(X,x)
%
%         X =   [n x p] data matrix for group of observations.
%         x =   [m x p] data matrix for observations, the distances of which are 
%                 to be determined.
%         --------------------------------------------------------------------
%         D2 =  [m x 1] vector of Mahalanobis distances.
%

% RE Strauss, 3/18/03

function D2 = mahaldist(X,x)
  if (isvector(x))
    x = x(:)';
  end;

  [n,p] = size(X);
  [m,q] = size(x);
  
  if (p~=q)
    error('  MahalDist: group and individual observations must have same dimensionality.');
  end;
  
  C = cov(X);
  centroid = mean(X);
  D2 = zeros(m,1);
  for i = 1:m
% i
% x_i = x(i,:)
    d = centroid-x(i,:);
    D2(i) = d * (C \ d');
  end;

  return;
  