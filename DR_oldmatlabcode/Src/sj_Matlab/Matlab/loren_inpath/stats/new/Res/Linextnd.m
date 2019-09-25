% LINEXTND: Extends a given line segment by a specified distance.
%
%     Usage: newpt = linextnd(pt1,pt2,dist)
%
%           pt1, pt2 = corresponding [n x 2] matrices x,y coordinates of n pairs of 
%                        points defining n lines, which are to be extended beyond pt2.
%           dist =     distance(s) for line to be extended; either an [n x 1] vector
%                        corresponding to pt1 & pt2, or a scalar to be applied to each
%                        of the n line segments.
%           -------------------------------------------------------------------------
%           newpt =    matching [n x 2] matrix of x,y coordinates of new endpoints of 
%                        n line segments.
%

function newpt = linextnd(pt1,pt2,dist)
  [n1,p1] = size(pt1);
  [n2,p2] = size(pt2);
  [nd,pd] = size(dist);

  if (p1~=2)                        % Check dimensions of input matrices
    if (p1==1 & n1==2)
      pt1 = pt1';
      [n1,p1] = size(pt1);
    else
      error('  Linextnd: input matrices of wrong dimension(s)');
    end;
  end;
  if (p2~=2)                        % Check dimensions of input matrices
    if (p2==1 & n2==2)
      pt2 = pt2';
      [n2,p2] = size(pt2);
    else
      error('  Linextnd: input matrices of wrong dimension(s)');
    end;
  end;

  if (n1>1 & n1~=n2)
    error('  Linextnd: dimensions of input matrices must agree');
  else
    n = n1;
  end;

  if (nd~=n)
    if (nd==1)
      dist = dist * ones(n,1);
    else
      error('  Linextnd: dimensions of input matrices must agree');
    end;
  end;

  newpt = zeros(n,2);

  for i = 1:n                     % Cycle thru point pairs
    dt = pt2(i,:) - pt1(i,:);
    d = sqrt(sum(dt.*dt));          % Distance between given points
    r = dist(i)./d;                 % Ratio of extension distance to that between pts
    newpt(i,:) = pt2(i,:) + dt*r;    % New point
  end;

  return;
