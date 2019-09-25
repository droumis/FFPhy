% PATHPLOT: Plots a path, with equal axes and nodes indicated.
%
%     Usage: pathplot(path,nonodes)
%
%           path =   [n x 2] matrix of point coordinates for a 2D path.
%           nondes = suppress plotting of path nodes.
%

% RE Strauss, 8/25/98
%   11/25/99 - changed plot colors for Matlab v5.

function pathplot(path,nonodes)
  if (nargin < 2) nonodes = []; end; 

  if (isempty(nonodes))
    nonodes = 0;
  end;

  x = path(:,1);
  y = path(:,2);

  plot(x,y,'k');
  if (~nonodes)
    hold on;
    plot(x,y,'ko');
    hold off;
  end;
  putbnd(x,y);
  axis('equal');
  
  return;

