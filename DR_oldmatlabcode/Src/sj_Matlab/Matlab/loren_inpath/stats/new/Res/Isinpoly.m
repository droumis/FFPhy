% ISINPOLY: Determines whether a point lies within a polygon.
%
%     Syntax: isin = isinpoly(pt,poly)
%
%           pt =   [n x 2] matrix of point coordinates.
%           poly = [m x 2] matrix specifying a closed polygon.
%           -----------------------------------------------------
%           isin = [n x 1] vector of boolean values: true (=1) if 
%                   corresponding point is inside the polygon, 
%                   false (=0) if outside.
%

% Sedgewick, R. 1988  Algorithms, 2nd ed. Addison Wesley (p. 353-355).

% RE Strauss, 5/4/97

function isin = isinpoly(pt,poly)
  X = 1; Y = 2;
  tol = 1e-6;

  [npts,p] = size(pt);
  if (p~=2)
    if (npts==2)
      pt = pt';
      [npts,p] = size(pt);
    else
      error('  Point matrix dimensions incorrect');
    end;
  end;

  [n,p] = size(poly);
  minpoly =  min(poly);
  maxpoly =  max(poly);
  meanpoly = mean(poly);

  if (poly(1,:)~=poly(n,:))   % Close polygon if open
    poly = [poly; poly(1,:)];
  else
    n = n-1;
  end;

  % List of line segments (x1,y1,x2,y2) comprising polygon

  polyedge = [poly(1:n,:) poly(2:(n+1),:)];

  % Find points lying within horizontal and vertical bounds of polygon

  inbounds = (pt(:,X) >= minpoly(X)-tol & pt(:,X) <= maxpoly(X)+tol & ...
              pt(:,Y) >= minpoly(Y)-tol & pt(:,Y) <= maxpoly(Y)+tol);

  isin = zeros(npts,1);       % Allocate output vector

  for ip = 1:npts             % Iterate thru points
    if (inbounds(ip))

      % Check if point is on boundary or vertex

      if (any(isonline(pt(ip,:),polyedge)))
        isin(ip)=1;
      else

        % Construct horizontal line segment from point to horizontal bound

        line = [pt(ip,X) pt(ip,Y) 0 pt(ip,Y)]; % Point is one end of line

        if (pt(ip,X) > meanpoly(X))         % If pt is to right of centroid,
          line(3) = minpoly(X)-1;          %   other end is before min x
        else
          line(3) = maxpoly(X)+1;          %   else other end is beyond max x
        end;

        % Jiggle line vertically toward centroid

        if (pt(ip,Y) > meanpoly(Y))
          line([2,4]) = line([2,4]) - (1e-12);
        else
          line([2,4]) = line([2,4]) + (1e-12);
        end;
        
        % Determine number of intersections with polygon
        % Point is in if #intersections is odd

        count = sum(intrsect(line,polyedge));
        if (rem(count,2)==1)
          isin(ip) = 1;                  
        end;
      end;
    end;
  end;

  return;

