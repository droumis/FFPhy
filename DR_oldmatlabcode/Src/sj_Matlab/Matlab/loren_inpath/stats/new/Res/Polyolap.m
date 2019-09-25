% POLYOLAP:  Finds the area of overlap (intersection) of two 2-dimensional 
%            polygons, which need not be convex.
%
%     Usage:  [abs_area,rel_area] = polyolap(P,Q)
%
%            P = [np x 2] matrix of vertex coordinates for polygon P.
%            Q = [nq x 2] matrix of vertex coordinates for polygon Q.
%            --------------------------------------------------------
%            abs_area = absolute area of overlap.
%            rel_area = relative area of overlap.
%                                                  

function [abs_area,rel_area] = polyolap(P,Q)
  [np,ca] = size(P);
  [nq,cb] = size(Q);
  
  if (ca~=2 | cb~=2)
    error('  Polygon coordinate matrices must each have two columns');
  end;

  get_rel = 0;
  if (nargout>1)
    get_rel = 1;
  end;

  if (any(P(1,:)~=P(np,:)))    % Close polygons if open
    P = [P;P(1,:)];
  else
    np = np-1;
  end;
  if (any(Q(1,:)~=Q(nq,:)))
    Q = [Q;Q(1,:)];
  else
    nq = nq-1;
  end;

  % Determine the vertices of each polygon that lie within the other, 
  % and the intersection points of the polygons.  
  % Make list of subpolygon vertex ordinates, sorted vertically.

  PinQ = [isinpoly(P(1:np,:),Q);0];
  QinP = [isinpoly(Q(1:nq,:),P);0];

  P_edges = [P(1:np,:), P(2:np+1,:)];
  Q_edges = [Q(1:nq,:), Q(2:nq+1,:)];

  [intsct,x,y] = intrsect(P_edges,Q_edges);
  ipt = find(intsct(:));
  x = x(:);
  y = y(:);

  Y = sort([P(PinQ,2); Q(QinP,2); y(ipt)]);

  % Find heights and vertical address (y value) of lateral 
  % midpoints of trapezoids.  Find all horizontal intersections 
  % of midpoints with both polygons.

  xmin = min([P(:,1);Q(:,1)]);      % Extent of abscissa
  xmax = max([P(:,1);Q(:,1)]);

  M = [];
  for i = 1:(length(Y)-1)
    h = Y(i+1)-Y(i);                % Trapezoid height
    if (h > (1e-6))                   % If new trapezoid, save lateral midpoints
      y = (Y(i)+Y(i+1))/2;
      [intsct,x,y] = intrsect([xmin y xmax y],[P_edges;Q_edges]);
      x = x(intsct)';
      y = y(intsct)';
      [x,j] = sort(x);
      xy = [x y(j)];                  % Midpoints
      lx = length(x);
      ident = [prod((xy(1:lx-1,:)==xy(2:lx,:))'),0];
      xy = xy(~ident,:);              % Delete repeated points
      M = [M; xy h*ones(size(xy,1),1)]; % Keep midpoints & corresponding height
    end;
  end;

  % Keep only the horizontal intersections that are common to both 
  % polygons.  These are the midpoints on the subpolygons.

  midpoint = M((isinpoly(M(:,1:2),P) & isinpoly(M(:,1:2),Q)),:);

  abs_area = 0;                           % Absolute overlap
  for i = 1:2:size(midpoint,1)
    w = midpoint(i+1,1)-midpoint(i,1);
    abs_area = abs_area + w*midpoint(i,3);
  end;

  if (get_rel)                            % Relative overlap
    rel_area = 2*abs_area/(polyarea(P)+polyarea(Q));
  end;

  return;

