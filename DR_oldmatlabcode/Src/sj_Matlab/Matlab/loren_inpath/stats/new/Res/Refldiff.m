% REFLDIFF: Given a polygon (specified by the point configuration) and the an axis 
%           specified by slope and intercept, finds the difference in outlines 
%           between the original polygon and its reflection about the axis.
%
%     Usage: diff = refldiff(axis,pts)
%
%           axis = [b1, b0] vector of coefficients.
%           pt =   [p x 2] matrix of point coordinates.
%           ---------------------------------------------------
%           diff = cumulative difference between half-outlines.
%

function diff = refldiff(axis,pts)
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

  bound1 = pts(ipt1:ipt2,:);              % Split polygon into two portions,
  bound2 = [];                            %   above and below x-axis
  for i = (ipt1+npts):-1:ipt2
    bound2 = [bound2; pts(wrap(i,npts),:)];
  end;

  bound1(:,2) = abs(bound1(:,2));         % Flip below-axis to above-axis
  bound2(:,2) = abs(bound2(:,2));

  n1 = size(bound1,1);                    % Numbers of points on the two boundaries
  n2 = size(bound2,1);

  seg1 = [bound1(1:n1-1,:) bound1(2:n1,:)];   % Line segments of boundaries
  seg2 = [bound2(1:n2-1,:) bound2(2:n2,:)];

  i = find(seg1(:,1) > seg1(:,3));        % Reverse endpoints for segments not 
  seg1(i,:) = [seg1(i,3:4) seg1(i,1:2)];  %   traversing left to right
  i = find(seg2(:,1) > seg2(:,3));
  seg2(i,:) = [seg2(i,3:4) seg2(i,1:2)];

  xmin = min([bound1(:,1); bound2(:,1)]); % Find extent of projection onto x-axis
  xmax = max([bound1(:,1); bound2(:,1)]);
  ymax = max([bound1(:,1); bound2(:,1)]); % Max height

  diff = 0;
  intervals = 50;

  for x = linspace(xmin,xmax,intervals)   % Divide into equal intervals
                                            % Locate segments at this x value
    i1 = find(seg1(:,1)<=x & seg1(:,3)>=x);     
    i2 = find(seg2(:,1)<=x & seg2(:,3)>=x);

    y1 = 0;
    y2 = 0;

    if (~isempty(i1))
      [i,xi,y1] = intrsect([x 0 x ymax],seg1(i1,:));    % Find intersections with boundaries
    end;
    if (~isempty(i2))
      [i,xi,y2] = intrsect([x 0 x ymax],seg2(i2,:));
    end;

    diff = diff + (sum(y1)-sum(y2)).^2;   % Sum the squared differences between boundaries
  end;

  return;

