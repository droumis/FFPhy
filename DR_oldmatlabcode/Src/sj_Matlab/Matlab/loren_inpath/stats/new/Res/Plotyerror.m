% PLOTYERROR: Produces scatterplot with vertical errorbars.
%
%     Usage: plotyerror(x,y,e,{cbnds},{connect})
%
%         x,y =     [n x 1] vectors of data values.
%         e =       [n x 1] corresponding vector of vertical deviations (if 
%                     symmetric), or [n x 2] matrix of lower and upper vertical 
%                     deviations (if asymmetric).
%         cbnds =   optional boolean flag indicating, if true, that 'e' 
%                     represents confidence bounds rather than deviations from 
%                     y [default = 0].
%         connect = optional boolean flag indicating, if true, that data 
%                     values are to be connected [default = 0].
%

% RE Strauss, 10/23/00
%   6/5/01 -   allow for asymmetric errorbars and for connecting data values;
%                sort input by x if connecting data values.
%   11/26/01 - replaced range(x) with max(x)-min(x).

function plotyerror(x,y,e,cbnds,connect)
  if (nargin < 4) cbnds = []; end;
  if (nargin < 5) connect = []; end;

  if (isempty(cbnds))
    cbnds = 0;
  end;
  if (isempty(connect))
    connect = 0;
  end;

  if (~isvector(x) | ~isvector(y))
    error('  PLOTYERROR: input matrices x & y must be vectors.');
  else
    x = x(:);
    y = y(:);
  end;

  n = length(x);
  if (isvector(e))
    if (cbnds)
      error('  PLOTYERROR: confidence bounds require lower and upper values');
    end;
    e = e(:)*ones(1,2);
  end;
  lene = size(e,1);
  if (length(y)~=n | lene~=n)
    error('  PLOTYERROR: input matrices not compatible.');
  end;

  if (cbnds)
    e = abs(e-y*ones(1,2));
  end;

  if (connect)
    [x,i] = sort(x);
    y = y(i);
    e = e(i,:);
  end;

  ticklen = 0.02*(max(x)-min(x))/2;

  plot(x,y,'ko');
  if (connect)
    hold on;
    plot(x,y,'k');
    hold off;
  end;
  hold on;
  for i = 1:n
    if (isvector(e))
      if (e(i)>0)
        plot([x(i) x(i)],[y(i)-e(i) y(i)+e(i)],'k');
        plot([x(i)-ticklen x(i)+ticklen],[y(i)-e(i) y(i)-e(i)],'k');
        plot([x(i)-ticklen x(i)+ticklen],[y(i)+e(i) y(i)+e(i)],'k');
      end;
    else
      if (e(i,1)>0)
        plot([x(i) x(i)],[y(i)-e(i,1) y(i)],'k');
        plot([x(i)-ticklen x(i)+ticklen],[y(i)-e(i,1) y(i)-e(i,1)],'k');
      end;
      if (e(i,2)>0)
        plot([x(i) x(i)],[y(i)+e(i,2) y(i)],'k');
        plot([x(i)-ticklen x(i)+ticklen],[y(i)+e(i,2) y(i)+e(i,2)],'k');
      end;
    end;
  end;
  hold off;
  putbnd([x; x],[y-e(:,1); y+e(:,2)]);

  if (n<=10)
    puttick(1:n);
  end;

  return;
