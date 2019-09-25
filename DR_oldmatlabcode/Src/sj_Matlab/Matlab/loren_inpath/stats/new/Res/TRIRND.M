% TRIRND: Generate a random triangulation of a closed polygon via recursion.
%         Passed to the procedure is the complete binary edge-use matrix,
%         and a vector indicating the indices for a subpolygon.
%         Returns the updated edge matrix.
%
%     Usage: edge = trirnd(edge,polygon)
%

% RE Strauss, 10/26/94

function edge = trirnd(edge,polygon)
  N = max(size(edge));                % Number of polygon points
  n = max(size(polygon));             % Number of subpolygon points

  % Randomly choose two points separated by at least one other point on
  % both sides.  Connect with an edge.
  i = 0;
  j = 0;
  while (any([(j-i<=1),(wrap(j+1,n)==i)]))
    R = floor(n*rand(1,2)+1);
    i = min(R);
    j = max(R);
  end;
  edge(polygon(i),polygon(j)) = 1;
  edge(polygon(j),polygon(i)) = 1;

  % Divide polygon into two parts, maintaining the original point labels.
  v = zeros(size(polygon));
  v(i:j) = ones(size(i:j));
  v1 = polygon(find(v));

  v = ones(size(polygon));
  v(i+1:j-1) = zeros(size(i+1:j-1));
  v2 = polygon(find(v));

  % Recursively triangulate each part that contains more than three points.
  if (max(size(v1)) > 3)
    edge = trirnd(edge,v1);
  end;
  if (max(size(v2)) > 3)
    edge = trirnd(edge,v2);
  end;

  return;




