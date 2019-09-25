% HOMERANGEF:  Objective function for homerange().
%
%     Usage: hr = homerangef(crds,nu1,nu2,nu3,stat)
%
%             crds =  [n x 2] matrix of coordinates of points of capture.
%             nu1,nu2,nu3 = not used (sent by bootstrp).
%             stat =  statistic used to estimate home range:
%                       1 = area of convex hull;
%                       2 = mean squared deviation from centroid;
%                       3 = mean squared deviation from 50% trimmed centroid;
%                       4 = median squared deviation from medioid.
%             ---------------------------------------------------------------
%             hr =    home-range estimated.
%                       

% RE Strauss, 7/26/00

function hr = homerangef(crds,nu1,nu2,nu3,stat)
  switch(stat)
    case 1,                             % Area of convex hull
      hr = polyarea(hull(crds));

    case 2,                             % MS deviation from centroid
      cent = mean(crds);
      hr = mean(eucl(crds,cent).^2);

    case 3,                             % MS deviation from trimmed centroid
      cent = centroid(crds);
      hr = mean(eucl(crds,cent).^2);

    case 4,                             % Median squared deviation from medoid
      [cent,med] = centroid(crds);
      hr = median(eucl(crds,med).^2);

    otherwise
      error('  HOMERANGEF: invalid statistic identifier.');
  end;

  return;
