% AreaInPoly1: Finds area for AreaInPoly()with 1 corner in polygon.
%
%     Usage: A = areainpoly1(cellaxes,poly)
%
%         cellaxes =  [2 x 4] matrix of pt coordinates of grid edges intersecting
%                       the polygon.
%         poly =      [n x 2] matrix of polygon vertices.
%         -----------------------------------------------------------------------
%         A =         area of grid cell lying within poloygon.
%

% RE Strauss, 5/20/03

function A = areainpoly1(cellaxes,poly)
  [intsct,xx,yy] = intrsect(cellaxes,poly);
  i1 = find(intsct(1,:));
  i1 = i1(1);
  i2 = find(intsct(2,:));
  i2 = i2(end);
  if (i1==i2)
    p = [cellaxes(1,1:2);xx(1,i1),yy(1,i1);xx(2,i2),yy(2,i2);cellaxes(1,1:2)];
  else
    if (i1<i2)
      s = i1+1:i2;
    else
      s = [(i1+1):(size(poly,1)-1),1:i2];
    end;
    p = [cellaxes(1,1:2);xx(1,i1),yy(1,i1);poly(s,:); ...
         xx(2,i2),yy(2,i2);cellaxes(1,1:2)];
  end;
  A = polyarea(p);
  
  return;
  
