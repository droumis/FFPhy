% LINEEQN: Given the coordinates of two points, calculates the intercept and 
%       slope of the line connecting them.
%
%     Usage: b = lineeqn(pt1,pt2)
%
%       pt1 = vector of x,y coordinates of first point.
%       pt2 = vector of x,y coordinates of second point.
%       ------------------------------------------------
%       b =   [2,1] vector of intercept and slope.
%

% RE Strauss, 1/20/98

function b = lineeqn(pt1,pt2)
  if (length(pt1)~=2 | length(pt2)~=2)
    error('  LINEEQN: input vectors of wrong length');
  end;    

  b = zeros(2,1);                  % Allocate return vector
  diff = pt2 - pt1;
  b(2) = diff(2)./diff(1);
  b(1) = pt2(2) - (pt2(1)*b(2));

  return;
