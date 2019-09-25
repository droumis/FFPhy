% REFLSTAT: Given a polygon (specified by the point configuration) and the an axis 
%           specified by slope and intercept, finds the differences in area and 
%           perimeter between the two sides.
%
%     Usage: [asym_area,asym_perim] = refldiff(axis,pts)
%
%           axis = [b1, b0] vector of coefficients.
%           pt =   [p x 2] matrix of point coordinates.
%           ---------------------------------------------------
%           asym_area = difference in half-areas.
%           asym_perim = difference in half-perimeters.
%

function [asym_area,asym_perim] = refldiff(axis,pts)
  [npts,c] = size(pts);
  b1 = axis(1);
  b0 = axis(2);

  if (pts(1,:) ~= pts(npts,:))            % Close polygon if necessary
    pts = [pts; pts(1,:)];
    npts = npts+1;
  end;

  xmin = min(pts(:,1));                   % Convert axis to line segment
  xmax = max(pts(:,1));                   %   extending thru polygon
  a1 = [2*xmin, 2*b1*xmin+b0];
  a2 = [2*xmax, 2*b1*xmax+b0];
  
  [intsct,xint,yint] = intrsect([a1 a2],[pts(1:npts-1,:) pts(2:npts,:)]);
  ipt = find(intsct);                     % Find intersections of axis with polygon

  if (length(ipt) < 1)
    diff = 1/eps;
    return;
  end;

  ipt1 = ipt(1);
  ipt2 = ipt(2);

  xint = xint(ipt)';                      % Intersection coordinates
  yint = yint(ipt)';

  p = [pts(1:ipt1,:); 
       xint(1) yint(1); 
       pts(ipt1+1:ipt2,:); 
       xint(2) yint(2);
       pts(ipt2+1:npts,:)];
  ipt1 = ipt1+1;
  ipt2 = ipt2+2;
  npts = npts+1;
  p(npts+1,:) = [];                       % Open polygon

  pts = regrot(ipt1,ipt2,p);              % Register and rotate polygon

  bound2 = pts(ipt1:ipt2,:);              % Split polygon into two portions,
  bound1 = [];                            %   above and below x-axis
  for i = (ipt1+npts):-1:ipt2
    bound1 = [bound1; pts(wrap(i,npts),:)];
  end;

  [area1,perim1] = polyarea(bound1);      % Stats for the two sides
  [area2,perim2] = polyarea(bound2);

  asym_area = area1 - area2;
  asym_perim = perim1 - perim2;

  return;

