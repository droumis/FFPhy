% GRPHULLS: Given a group-identification vector and matrix of coordinates,
%           returns a matrix of hulls delimited by row vectors of zeros,
%           and a matrix of group centroids.
%
%     Usage: [hulls, centroids] = grphulls(grps,crds)
%
%         grps =      [N x 1] group-identification vector.
%         crds =      [N x 2] matrix of X,Y coordinates.
%         ---------------------------------------------------------------
%         hulls =     [h x 2] matrix of hulls delimited by rows of zeros.
%         centroids = [g x 2] matrix of group centroids.
%

% RE Strauss, 8/9/94

function [hulls,centroids] = grphulls(grps,crds)
  G = design(grps);                     % Design matrix
  [N,ngrps] = size(G);
  hulls = [];
  centroids = [];

  for g = 1:ngrps
    index = find(G(:,g));               % Indices to nonzero elements
    points = crds(index,:);             % Points for current group
    h = hull(points);                   % Hull for current group
    hulls = [hulls; h; 0 0];            % Accumulate hulls
    c = mean(points);                   % Group centroid
    centroids = [centroids; c];         % Accumulate centroids
  end;

  [H,c] = size(hulls);
  hulls(H,:) = [];                      % Delete final delimiter

  return;
