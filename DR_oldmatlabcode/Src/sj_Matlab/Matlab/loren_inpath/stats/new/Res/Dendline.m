% DENDLINE: Indicates groups of taxa on an existing UPGMA dendrogram by drawing 
%             vertical lines.  Assumes that the UPGMA is in the current active figure.
%
%     Usage: dendline(x,bounds)
%
%           x =       number of units, to the right of the zero position on the X axis,
%                       at which lines are to be drawn.
%           bounds = [n x 2] matrix specifying the taxon-position bounds for which lines 
%                       are to be drawn.  Taxa are number consecutively from bottom to 
%                       top of the dendrogram.
%

% RE Strauss, 6/8/98
%   11/25/99 - changed plot colors for Matlab v5.

function dendline(x,bounds)
  [n,p] = size(bounds);
  x = -x;

  ext = 0.25;

  v = axis;                 % Extends the X axis
  v(1) = x;
  axis(v);

  hold on;
  plot([eps x],[0 0],'k');
  plot([eps x],[v(4) v(4)],'k');

  for i = 1:n
    y1 = bounds(i,1) - ext;
    y2 = bounds(i,2) + ext;
    plot([x x],[y1 y2],'k');
  end;
  
  hold off;
  return;
