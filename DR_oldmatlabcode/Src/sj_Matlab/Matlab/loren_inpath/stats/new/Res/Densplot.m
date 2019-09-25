% DENSPLOT: Scatterplot of the densities of two discrete variables, 
%           wherein each point on the grid indicates the number of 
%           observations.
%
%     Syntax:  densplot(x,y,{fontsize})
%

% RE Strauss, 3/9/97
%   9/19/99 - updated handling of default input arguments.

function densplot(x,y,fontsize)
  if (nargin<3) fontsize = []; end;

  if (min(size(x))>1)
    error('  Input not in vector form.');
  end;
  if (size(x) ~= size(y))
    error('  Input vectors not compatible.');
  end;

  n = length(x);
  if (isempty(fontsize))
    fontsize = 10;
  end;

  xmin = min(x)-1;
  ymin = min(y)-1;
  dens = zeros(range(y)+1,range(x)+1);

  for i = 1:n
    c = x(i)-xmin;
    r = y(i)-ymin;
    dens(r,c) = dens(r,c)+1;
  end;

  plot(x,y,'k.');
  putbnd(x,y);

  [r,c] = size(dens);
  dx = xmin-0.01*range(x);
  dy = ymin;

  hold on;
  for i = 1:c
    for j = 1:r
      d = dens(j,i);
      if (d>0)
        h = text(i+dx,j+dy,num2str(d));
%        get(h)
        set(h,'FontSize',fontsize,'Color',[1 1 1])
      end;
    end;
  end;

  hold off;

  return;
