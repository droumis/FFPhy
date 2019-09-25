% BOOTMEAN: Given a bivariate distribution, plots the centroid and convex 
%           hull and and a distribution of boostrapped estimates of the mean.
%
%       Usage: bootmean(X,nboot)
%
%           X -     [n x 2] matrix of data points
%           nboot - number of bootstrapped means to be plotted
%

% RE Strauss, 7/31/95
%   8/19/99 - changed plot colors for Matlab v5

function bootmean(X,nboot)
  c = mean(X)                          % Centroid
  hull_pts = hull(X);                   % Convex hull
  csave = zeros(nboot,2);               % Save bootstrapped means

  hold on;
  plot(hull_pts(:,1),hull_pts(:,2),'k');    % Plot hull
  axis('square');
  plot(c(1),c(2),'*b');

  for i=1:nboot
    Y = bootsamp(X);
    c = mean(Y);
    plot(c(1),c(2),'+r');
    csave(i,:) = c;
  end;
  grandc = mean(csave)

  hold off;
  
  return;
