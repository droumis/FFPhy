% MINSIZE: Select data for groups having at least a minimum sample size.
%
%     Usage: [rX,rg] = minsize(X,grps,minsampsize)
%
%         X =       [n x p] data matrix.
%         grps =    [n x 1] group-membership vector.
%         minsize = minimum sample size for inclusion.
%         --------------------------------------------
%         rX =      reduced [m x p] data matrix.
%         rg =      reduced [m x 1] grouping vector.
%

% RE Strauss, 8/3/00

function [rX,rg] = minsize(X,grps,minsampsize)
  [n,p] = size(X);
  ng = length(grps);

  if (ng ~= n)
    error('  MINSIZE: data matrix and grouping vector not of same size.');
  end;

  [u,f] = uniquef(grps);
  incl = find(f>=minsampsize);

  rX = [];
  rg = [];

  for i = 1:length(incl)
    g = u(incl(i));
    j = find(grps == g);
    rX = [rX; X(j,:)];
    rg = [rg; grps(j)];
  end;

  return;
