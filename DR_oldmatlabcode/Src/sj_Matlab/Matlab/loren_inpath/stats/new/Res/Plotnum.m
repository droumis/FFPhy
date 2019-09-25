% Plotnum: Modification of plot() command to plot observation numbers 
%          rather than symbols.  If one input matrix is passed rather than
%          two vectors, the first two columns of the matrix are plotted.
%
%     Syntax:  plotnum(x,y,{fontsize})
%                 OR
%              plotnum([x,y],{fontsize})
%

% RE Strauss, 9/6/96
%   1/2/99 -  make second argument optional.
%   9/7/99 -  changed plot colors for Matlab v5.
%   3/21/00 - simplify plot and buffer logic.
%   10/3/00 - produce initial plot with white points so that plot can 
%               subsequentaly be rescaled with axis().
%   6/1/01 -  plot points as black dots on figure 1.

function plotnum(x,y,fontsize)
  if (nargin<2) y = []; end;
  if (nargin<3) fontsize = []; end;

  if (isempty(y) | isscalar(y))
    fontsize = y;
    if (~isvector(x))
      y = x(:,2);
      x = x(:,1);
    else
      error('  PLOTNUM: invalid coordinate vectors');
    end;
  end;

  if (isempty(fontsize))
    fontsize = 10;
  end;

  if (length(x)~=length(y))
    error('  PLOTNUM: coordinate vectors must be of same length');
  end;

  deltax = 0.011 * range(x);

  plot(x,y,'k.');
  hold on;
  for i = 1:length(x)
    h = text(x(i)+deltax,y(i),int2str(i));
    set(h,'fontsize',fontsize);
  end;
  hold off;
  putbnd([x,y],0.06);

  return;
