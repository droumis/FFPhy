% LineCrds: Given the two endpoints of a line segment, or the intercept and slope of a line,
%           returns the line coordinates [l1,l2,l3] where l1p1 + l2p2 + l3 = 0 for any point
%           [p1,p2] on the line.
%
%     Usage: line = linecrds(p,q)    OR
%            line = linecrds([p,q])  OR
%            line = linecrds(b0,b1)  OR
%            line = linecrds([b0,b1])
%
%         p,q = 2-element vectors of point coordinates.
%         b0,b1 = scalar values for intercept and slope.
%         ----------------------------------------------
%         line = 3-element vector of line coordinates.
%

% Bookstein et al. 1985. Morphometrics in Evolutionary Biology, appendix A.1.

% RE Strauss, 4/22/03

function line = linecrds(a,b)
  a = a(:);
  lena = length(a);
  if (nargin == 2)
    b = b(:);
    lenb = length(b);
  end;
  
  
  err = 0;
  if (lena~=1 & lena~=2 & lena~=4)
    err = 1;
  end;
  if (nargin == 2)
    if (lenb~=1 & lenb~=2 & lenb~=4)
      err = 1;
    end;
  end;
  if (err)
    error('  LineCrds: invalid input arguments');
  end;

  pq = 0;
  
  if (nargin == 1)
    if (lena==2)
      b0 = a(1);
      b1 = a(2);
    else
      p = a(1:2);
      q = a(3:4);
      pq = 1;
    end;
  else %if (nargin == 2)
    if (lena==2)
      p = a;
      q = b;
      pq = 1;
    else
      b0 = a;
      b1 = b;
    end;
  end;
  
  if (pq)
    if (p(1)==q(1))               % Vertical line
      line = [-1,0,1];
    else
      b1 = (p(2)-q(2))/(p(1)-q(1));
      b0 = p(2) - b1*p(1);
      line = [b1,-1,b0];
    end;
  else
    line = [b1,-1,b0];
  end;

  return;
  
