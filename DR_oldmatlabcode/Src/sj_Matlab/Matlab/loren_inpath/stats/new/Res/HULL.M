% HULL: From a set of point coordinates, find convex hull.  Returns only hull
%       vertex points, disregarding non-vertex hull points.
%
%     Usage: [hull_pts,index] = hull(crds)
%
%        crds =       [n x 2] matrix of x,y point coordinates
%        ---------------------------------------------------------------
%        hull_pts =   [h x 2] matrix of x,y hull coordinates
%        hull_index = [h x 1] vector of indices to hull points in 'crds'
%

% Sedgewick, R. 1988  Algorithms, 2nd ed. Addison Wesley (Ch. 25).
%   Uses the package-wrapping algorithm.

% RE Strauss, 1/5/96
%   10/4/01 - changed 'continue' variable due to v6 usage.

function [hull_pts,hull_index] = hull(crds)
  TRUE = 1; FALSE = 0;

  [n,p] = size(crds);
  if (p~=2)
    error('  HULL: only 2-dimensional hulls allowed');
  end;

  if (n<3)                            % Fewer than 3 points passed
    if (n==1)                           % Single point passed
      hull_index = 1;
      hull_pts = crds;
    else                                % Two points passed
      [m,miny] = min(crds(:,2));
      if (miny==1)
        hull_index = [1 2 1]';
      else
        hull_index = [2 1 2]';
      end;
      hull_pts = [crds(hull_index,:)];
    end;
    return;
  end;

  m = min(crds(:,2));                 % Find min y value ("canonical form")
  minpt = find(crds(:,2)==m);
  if (length(minpt)>1)                % If ties,
    xi = 1;                           %   find min x among set of min y's
    xival = crds(minpt(1),1);
    for i=2:length(minpt)
      if (crds(minpt(i),1) < xival)
        xi = i;
        xival = crds(minpt(i),1);
      end;
    end;
    minpt = minpt(xi);
  end;

  crds = [crds; crds(minpt,:)];       % Append min-y pt to end of coords
  seq = [1:n+1]';                     % Input sequence of points

  minangle = 0;
  maxangle = theta([0 0],[1 -eps]);
  nh = 0;

  contin = 1;
  while (contin)
    if (nh>1 & crds(minpt,:)==crds(1,:))
      contin = 0;
    else
      nh = nh+1;                      % Update current hull point
      tpt = crds(nh,:);               % Switch current and min-angle points
      tseq = seq(nh);
      crds(nh,:) = crds(minpt,:);
      seq(nh) = seq(minpt);
      crds(minpt,:) = tpt;
      seq(minpt) = tseq;

      [angle,ident] = theta(crds(nh,:),crds((nh+1):(n+1),:)); % Get angles from hull pt
      angle = angle - minangle;                               % With respect to prev angle

      while (min(angle)<0)                          % Wrap negative angles around
        indx = find(angle<0);                       %   to positive
        angle(indx) = angle(indx) + maxangle;
      end;

      indx = find(ident>0);                         % Disregard pts identical to hull pt
      angle(indx) = maxangle * ones(length(indx),1);

      [ma,minpt] = min(angle);                      % Find min angle and point
      prev_minangle = minangle;
      minangle = minangle + ma;
      minpt = minpt + nh;
    end;
  end;

  hull_index = [seq(1:nh); seq(1)];     % Stash point indices
  hull_pts = [crds(1:nh,:); crds(1,:)]; % Coordinates of hull points

%  i = 0;
%  while (i < nh-2)
%    i = i+1;
%    hp = hull_pts(i:(i+2),:);           % Check for colinear triplets of hull points
%    if (ptdir(hp)==0)                   % If found, delete the center one
%      hull_pts(i+1,:) = [];
%      hull_index(i+1) = [];
%      nh = nh-1;
%      i = i-1;
%    end;
%  end;

  return;
