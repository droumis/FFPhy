% PATHLEN:  Given a set of point coordinates defining a path (in any number of 
%           dimensions), finds the path length.  The path length of a single 
%           point or of a null matrix is defined to be zero.
%
%     Usage: len = pathlen(crds)
%
%           crds = [n x p] matrix of point coordinates.
%           -------------------------------------------
%           len =  total path length.
%

% RE Strauss, 9/14/99

function len = pathlen(crds)
  [n,p] = size(crds);

  if (n<2 & p>1)
    crds = crds';
    [n,p] = size(crds);
  end;

  if (n<2)
    len = 0;
  elseif (n==2)
    d = crds(2,:)-crds(1,:);
    if (p>1)
      len = sqrt(sum(d'.^2));
    else
      len = abs(crds(1)-crds(2));
    end;
  else
    c1 = crds(1:n-1,:);
    c2 = crds(2:n,:);
    d = c2-c1;

    if (p>1)
      len = sum(sqrt(sum(d'.^2))');
    else
      len = sum(d);
    end;
  end;

  return;
