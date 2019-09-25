% INTRSECT: Determines whether and where line segments intersect.  
%           Given two sets of line segments, the function does 
%           pairwise comparisons of all segments in the first set 
%           with all segments in the second set.  Optionally returns 
%           the X and Y coordinates of intersection points in 
%           separate matrices.  
%           
%     Usage: [intsct,x,y] = intrsect(line1,line2)
%
%           line1 =  [n1 x 4] matrix specifying the point-coordinates (x,y) 
%                      of n line-segment endpoints: [x1 y1, x2 y2].
%                         OR
%                    [n1+1 x 2] matrix of point-coordinates specifying a path of
%                      connected line segments.
%           line2 =  [n2 x 4] matrix specifying the corresponding 
%                      point-coordinates of m line-segment endpoints.
%                         OR
%                    [n2+1 x 2] matrix of point-coordinates specifying a path of
%                      connected line segments.
%           -----------------------------------------------------------
%           intsct = [n1 x n2] matrix of boolean variables, true (=1) 
%                      if segments intersect, false (=0) if not or if line 
%                      segments are coincident.
%           x =      [n1 x n2] matrix of corresponding x coordinates, with 
%                      zeros for non-intersecting pairs of segments.
%           y =      [n1 x n2] matrix of correspoinding y coordinates.
%

% Sedgewick, R. 1988  Algorithms, 2nd ed. Addison Wesley (p. 351).
% Bookstein et al. 1985. Morphometrics in Evolutionary Biology (p. 215-216).

% RE Strauss, 1/22/96
%   3/14/03 - allow alternate input of paths.

function [intsct,x,y] = intrsect(line1,line2)
  get_crd = 0;
  if (nargout > 1) get_crd = 1; end;

  [n1,p1] = size(line1);
  [n2,p2] = size(line2);

  if ((p1~=4 & p1~=2)| (p2~=4 & p2~=2))
    error('INTRSECT: Input matrices must have 2 or 4 columns');
  end;
  
  if (p1==2)
    line1 = [line1(1:end-1,:),line1(2:end,:)];
    n1 = n1-1;
  end;
  if (p2==2)
    line2 = [line2(1:end-1,:),line2(2:end,:)];
    n2 = n2-1;
  end;

  if (get_crd)
    x = zeros(n1,n2);
    y = x;
  end;

  intsct = zeros(n1,n2);        % Allocate boolean result matrix
  p1 = [1,2];                   % Coordinates of segment endpoints
  p2 = [3,4];

  for i = 1:n1                  % Cycle thru pairs of line segments
    l11 = line1(i,p1);
    l12 = line1(i,p2);

    for j = 1:n2
      l21 = line2(j,p1);
      l22 = line2(j,p2);

      dir1 = ptdir(l11,l12,l21);
      dir2 = ptdir(l11,l12,l22);
      dir3 = ptdir(l21,l22,l11);
      dir4 = ptdir(l21,l22,l12);
      intsct(i,j) = ((dir1*dir2 <= 0) & (dir3*dir4 <= 0));

      if (intsct(i,j) & get_crd)  % Get intersection coordinates
        % Reparameterize lines in terms of endpoints

        l1 = [l11(2)-l12(2), l12(1)-l11(1), l11(1)*l12(2)-l12(1)*l11(2)];
        l2 = [l21(2)-l22(2), l22(1)-l21(1), l21(1)*l22(2)-l22(1)*l21(2)];
        d = l1(1)*l2(2) - l1(2)*l2(1);

        if (abs(d) < eps)   % If d=0, lines are parallel or coincident
          intsct(i,j) = 0;
        else
          x(i,j) = (l1(2)*l2(3) - l1(3)*l2(2))/d;
          y(i,j) = (l1(3)*l2(1) - l1(1)*l2(3))/d;
        end;
      end;
    end;
  end;
  
  return;

