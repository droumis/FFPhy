% Arrow: Plots an arrow from one specified point to another on the current plot.  
%        Axis boundaries and aspect ratio of plot should be set before arrow() 
%        is called.
%
%     Syntax:   arrow(tail,head,{issquare},{headsize},{style})
%
%           tail =      vector of x,y coordinates of tail point.
%           head =      vector of x,y coordinates of head point.
%           issquare =  boolean flag indicating that the plot has been squared 
%                         [default = 0].
%           headsize =  size of arrow head (proportion of arrow shaft length) 
%                         [default = 5%].
%           style =     string vector indicating color/style of line 
%                         [default = black solid line].
%

% RE Strauss, 2/7/97
%   9/3/99 -    misc changes for Matlab v5.
%   11/25/99 -  changed default arrow-head size.
%   2/26/00 -   changed names of input arguments.

function arrow(tail,head,issquare,headsize,style)
  if (nargin < 3) issquare = []; end;
  if (nargin < 4) headsize = []; end;
  if (nargin < 5) style = []; end;

  headangle = 0.90*pi;
  
  stail = ptscale(tail);        % Scale coords to proportional distances
  shead = ptscale(head);

  if (isempty(headsize))
    headsize = 0.05;
  end;
  if (isempty(style))
    style = 'k-';
  end;
  if (isempty(issquare))
    issquare = 0;
  end;

  if (headsize<0 | headsize>1)
    error('  ARROW: arrowhead size must be expressed as proportion');
  end;

  if (issquare)
    aspect_ratio = 1;
  else
    rect = get(gca,'Position');
    aspect_ratio = (rect(4)-rect(2))./(rect(3)-rect(1));
  end;

  theta = angle((shead(1)-stail(1)) + (shead(2)-stail(2))*sqrt(-1));  % Shaft angle
  theta1 = theta + headangle;
  theta2 = theta - headangle;
  
  x1 = shead(1) + headsize * cos(theta1) / aspect_ratio;
  y1 = shead(2) + headsize * sin(theta1);
  p =  ptscale([x1,y1],1);
  x1 = p(1);
  y1 = p(2);

  x2 = shead(1) + headsize * cos(theta2) / aspect_ratio;
  y2 = shead(2) + headsize * sin(theta2);
  p =  ptscale([x2,y2],1);
  x2 = p(1);
  y2 = p(2);

  hold on;
  plot([tail(1) head(1)],[tail(2) head(2)],style);
  plot([head(1) x1],[head(2) y1],style);
  plot([head(1) x2],[head(2) y2],style);

  return;
