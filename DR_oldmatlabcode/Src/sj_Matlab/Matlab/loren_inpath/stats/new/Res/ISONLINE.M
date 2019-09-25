% ISONLINE: Determines whether a point lies on a line or line segment, 
%           within a specified tolerance.  The line (segment) is specified 
%           by two points.
%
%     Syntax: ison = isonline(pt,line,{extended},{tol})
%
%           pt =   [n x 2] matrix of coordinates for n points.
%           line = [m x 4] matrix specifying the endpoints (x1,y1,x2,y2) 
%                    for m line segments.
%           extended = boolean flag indicating that the point is to be 
%                    compared to the line passing thru (x1,y1) and (x2,y2) 
%                    rather than the line segment [default=0].
%           tol =  tolerance around line [default = 1e-6].
%           -----------------------------------------------------
%           ison = [n x m] matrix of boolean values: true (=1) if 
%                   the n_th point lies on the m_th line, 
%                   false (=0) if not.
%

function ison = isonline(pt,line,extended,tol)
  if (nargin < 3)
    extended = 0;
  end;
  if (nargin < 4)
    tol = 1e-6;
  end;

  [n,cn] = size(pt);
  [m,cm] = size(line);

  if (cn~=2 | cm~=4)
    error('  Input matrix is wrong size');
  end;

  ison = zeros(n,m);      % Allocate output matrix

  for i = 1:n             % Cycle thru points
    px = pt(i,1);
    py = pt(i,2);

    for j = 1:m             % Cycle thru lines
      p = line(j,1:2);
      q = line(j,3:4);

      l1 = p(2)-q(2);
      l2 = q(1)-p(1);
      l3 = p(1)*q(2)-p(2)*q(1);
      d = abs((px*l1+py*l2+l3)/(sqrt(l1*l1+l2*l2)));

      if (d <= tol)
        if (extended)
          ison(i,j) = 1;
        else
          if (sum(eucl(pt(i,:),[p;q]))-eucl(p,q) < 10*tol)
            ison(i,j) = 1;
          end;
        end;
      end;
    end;
  end;

  return;
