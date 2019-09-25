% PTDIR: Determines whether the path between three points p1,p2,p3 is
%        counterclockwise or clockwise.
%
%     Usage: path = ptdir(p1,p2,p3)
%
%        p1,p2,p3 = Matching [n x 2] matrices containing point coordinates.
%                     If only p1 is passed, it is assumed to be a [3 x 2]
%                     matrix specifying [p1;p2;p3].
%        --------------------------------------------------------------------
%        path =     [n x 1] vector of path directions: 
%                     1 if counterclockwise, -1 if clockwise, 0 if colinear.
%

% Sedgewick, R. 1988  Algorithms, 2nd ed. Addison Wesley (p. 350).
%   Modified so that colinear points always return a 0, regardless of their 
%   sequence.

% RE Strauss, 3/25/98

function path = ptdir(p1,p2,p3)
  if (nargin < 2)
    if (size(p1,1)==3)
      p2 = p1(2,:);
      p3 = p1(3,:);
      p1 = p1(1,:);
    else
      error('PTDIR: Invalid input arguments');
    end;
  end;

  [r1,c1] = size(p1);
  [r2,c2] = size(p2);
  [r3,c3] = size(p3);

  if (sum([c1 c2 c3])~=6)
    error('PTDIR: Invalid input arguments');
  end;

  if (r1~=r3 | r2~=r3)
    error('PTDIR: input matrices not compatible');
  end;

  path = zeros(r1,1);

  for i = 1:r1                        % For each set of points
    dx1 = p2(i,1) - p1(i,1);
    dy1 = p2(i,2) - p1(i,2);
    dx2 = p3(i,1) - p1(i,1);
    dy2 = p3(i,2) - p1(i,2);

    cond1 = dx1*dy2;
    cond2 = dy1*dx2;

    if (abs(cond1-cond2)<eps)           % Colinear
      path(i) = 0;
    else
      if (cond1 > cond2)                % Counterclockwise
        path(i) = 1;
      else                              % Clockwise
        path(i) = -1;
      end;
    end;
  end;

  return;

