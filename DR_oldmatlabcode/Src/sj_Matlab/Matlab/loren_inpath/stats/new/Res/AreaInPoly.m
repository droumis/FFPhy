% AreaInPoly: Given a set of points describing a closed polygon and the dimensions of a 
%             cartesian gridon which the polygon occurs, returns a matrix specifying the 
%             proportion of area of each grid cell that lies within the polygon.  
%             The polygon is assumed to be specified in counterclockwise sequence, and
%             is assumed to enter each grid cell only once.
%             Note that the grid position (0,0) is at lower left and (r,c) is at upper 
%             right, so that the resulting matrix is "coordinate mode": row 1 is at 
%             bottom, and row r is at top.
%             
%     Usage: G = areainpoly(poly,[xmax,ymax],{doplot})
%
%         poly =    [n x 2] matrix of polygon vertices.
%         dim =     vector (length 2) of numbers of columns (ceil(xmax)) and 
%                     rows (ceil(ymax))of matrix (grid) on which the polygon lies.
%         doplot =  optional boolean value indicating, if true, that a plot of grid,
%                     polygon, and cell values is to be produced [default = 0].
%         --------------------------------------------------------------------------
%         G =       matrix indicating proportion of area of each grid cell lying 
%                     within polygon.
%

% RE Strauss, 5/20/03

function G = areainpoly(poly,dim,doplot)
  if (~nargin) help areainpoly; return; end;
  
  if (nargin < 3) doplot = []; end;
  if (isempty(doplot)) doplot = 0; end;
  
  [n,p] = size(poly);
  if (p~=2)
    error('  AreaInPoly: polygon matrix must be size [n x 2].');
  end;
  
  c = ceil(dim(1));
  r = ceil(dim(2));
  G = zeros(r,c);
  
%   poly = makepolygon(poly);
  if (poly(1,:)~=poly(end,:))   % Close polygon if open
    poly = [poly; poly(1,:)];
  end;

  grid = [makegrps(0:c,r+1), makerepeatseqs(r+1,c+1)-1];  % Matrix of grid intersections
  isin = isinpoly(grid,poly);

  for x = 1:c                   % For each grid cell,
    for y = 1:r
      p1 = [x-1,y-1];             % Lower left corner of grid cell
      p2 = [x-1,y];               % Upper left
      p3 = [x,y-1];               % Lower right
      p4 = [x,y];                 % Upper right
      
      i = [find(grid(:,1)==(x-1) & grid(:,2)==(y-1)), ...
           find(grid(:,1)==(x-1) & grid(:,2)==y), ...
           find(grid(:,1)==x & grid(:,2)==(y-1)), ...
           find(grid(:,1)==x & grid(:,2)==y)];
      v = tostr(isin(i))';    % Boolean string, 0=vertex out, 1=vertex in
      switch (v)              % v = [lower left, upper left, lower right, upper right]
        case ['0000'],
          G(y,x) = 0;
        case ['0001'],
          cellaxes = [p4,p2; p4,p3];
          G(y,x) = AreaInPoly1(cellaxes,poly);
        case ['0010'],
          cellaxes = [p3,p4 ;p3,p1];
          G(y,x) = AreaInPoly1(cellaxes,poly);            
        case ['0011'],
          cellaxes = [p4,p2; p3,p1];
          G(y,x) = AreaInPoly2(cellaxes,poly);
        case ['0100'],
          cellaxes = [p2,p1 ;p2,p4];
          G(y,x) = AreaInPoly1(cellaxes,poly);            
        case ['0101'],
          cellaxes = [p2,p1; p4,p3];
          G(y,x) = AreaInPoly2(cellaxes,poly);
        case ['0110'],
          error('  AreaInPoly: currently illegal.');
        case ['0111'],
          cellaxes = [p1,p2;p1,p3];
          G(y,x) = 1-AreaInPoly1(cellaxes,poly);
        case ['1000'],
          cellaxes = [p1,p3; p1,p2];
          G(y,x) = AreaInPoly1(cellaxes,poly);
        case ['1001'],
          error('  AreaInPoly: currently illegal.');
        case ['1010'],
          cellaxes = [p3,p4; p1,p2];
          G(y,x) = AreaInPoly2(cellaxes,poly);
        case ['1011'],
          cellaxes = [p2,p4; p2,p1];
          G(y,x) = 1-AreaInPoly1(cellaxes,poly);
        case ['1100'],
          cellaxes = [p1,p3; p2,p4];
          G(y,x) = AreaInPoly2(cellaxes,poly);
        case ['1101'],
          cellaxes = [p3,p1; p3,p4];
          G(y,x) = 1-AreaInPoly1(cellaxes,poly);
        case ['1110'],
          cellaxes = [p4,p3; p4,p2];
          G(y,x) = 1-AreaInPoly1(cellaxes,poly);
        case ['1111'],
          G(y,x) = 1;
      end;
    end;
  end;
  
  if (doplot)
    figure;
    hold on;
    for i = 0:c
      plot([i i],[0,r],':k');
    end;
    for i = 0:r
      plot([0,c],[i i],':k');
    end;
    plot(poly(:,1),poly(:,2),'k',poly(:,1),poly(:,2),'ok');
    hold off;
    axis([0 c 0 r]);
    axis equal;
    axis off;
  end;
  
  [i,j] = find(G<0 | G>1);          % Temporary fix for paths that enter grid cells > once
  if (~isempty(i))
    for k = 1:length(i)
      G(i(k),j(k)) = 0.5;
    end;
  end;
  
  G = flipud(G);
  
  return;
  