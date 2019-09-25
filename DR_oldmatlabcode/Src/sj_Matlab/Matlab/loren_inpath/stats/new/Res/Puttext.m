% PUTTEXT: Adds text to a figure, with default font characteristics.  
%          Differs from text() in that the position of the text is given 
%          in proportions of axis length rather than fixed coordinates.
%
%     Syntax: puttext(x,y,textstring,axisflag,fontsize)
%
%           x,y -        proportional positions of x,y coordinates, in relation 
%                         to current axis ranges.
%           textstring - character string, in single quotes.
%           axisflag -   flag indicating scaling of axes:
%                           0 = X linear, Y linear  [default]
%                           1 = X log10,  Y linear  (semilogx)
%                           2 = X linear, Y log10   (semilogy)
%                           3 = X log10,  X log10   (loglog)
%           fontsize -   font size [default = 14]
%

% RE Strauss, 11/7/97

function puttext(x,y,textstring,axisflag,fontsize)
  if (nargin < 4) axisflag = []; end;
  if (nargin < 5) fontsize = []; end;

  if (isempty(axisflag))
    axisflag = 0;
  end;
  if (isempty(fontsize))
    fontsize = 14;
  end;

  v = axis;                         % Get current axes min,max

  xmin = v(1); 
  xmax = v(2);
  if (axisflag == 0 | axisflag == 2)
    x = xmin + x*(xmax-xmin);
  else
    xmin = log10(xmin);
    xmax = log10(xmax);
    x = 10^(xmin + x*(xmax-xmin));
  end;

  ymin = v(3); 
  ymax = v(4);
  if (axisflag == 0 | axisflag == 1)
    y = ymin + y*(ymax-ymin);
  else
    ymin = log10(ymin);
    ymax = log10(ymax);
    y = 10^(ymin + y*(ymax-ymin));
  end;

  h = text(x,y,textstring);
  set(h,'FontName','helvetica','FontSize',fontsize);

  return;

