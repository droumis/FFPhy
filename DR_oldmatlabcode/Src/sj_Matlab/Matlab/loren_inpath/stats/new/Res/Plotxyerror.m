% PLOTXYERROR: Produces scatterplot with horizontal and vertical errorbars.
%
%     Usage: plotxyerror(x,y,xe,ye,{cbnds})
%
%         x,y =     [n x 1] vectors of data values.
%         xe =      [n x 1] corresponding vector of horizontal deviations (if 
%                     symmetric), or [n x 2] matrix of lower and upper horizontal 
%                     deviations (if asymmetric).
%         ye =      [n x 1] corresponding vector of vertical deviations (if 
%                     symmetric), or [n x 2] matrix of lower and upper vertical 
%                     deviations (if asymmetric).
%         cbnds =   optional boolean flag indicating, if true, that 'xe' and 'ye'
%                     represent confidence bounds rather than deviations from 
%                     y [default = 0].
%

% RE Strauss, 6/17/01, modified from plotyerror()

function plotxyerror(x,y,xe,ye,cbnds)
  if (nargin < 5) cbnds = []; end;

  if (isempty(cbnds))
    cbnds = 0;
  end;

  if (~isvector(x) | ~isvector(y))
    error('  PLOTXYERROR: input matrices x & y must be vectors.');
  else
    x = x(:);
    y = y(:);
  end;
  if (size(xe) ~= size(ye))
    error('  PLOTXYERROR: input error matrices not compatible.');
  end;

  n = length(x);
  if (isvector(xe) | isvector(ye))
    if (cbnds)
      error('  PLOTXYERROR: confidence bounds require lower and upper values');
    end;
    xe = xe(:)*ones(1,2);
    ye = ye(:)*ones(1,2);
  end;
  lene = size(xe,1);
  if (length(y)~=n | lene~=n)
    error('  PLOTXYERROR: input matrices not compatible.');
  end;

  if (cbnds)
    xe = abs(xe-y*ones(1,2));
    ye = abs(ye-y*ones(1,2));
  end;

  xticklen = 0.02*range(y)/2;
  yticklen = 0.02*range(x)/2;

  plot(x,y,'ko');
  hold on;
  for i = 1:n
    if (isvector(xe))
      if (xe(i)>0)
        plot([x(i)-xe(i) x(i)+xe(i)],[y(i) y(i)],'k');
        plot([x(i)-xe(i) x(i)-xe(i)],[y(i)-xticklen y(i)+xticklen],'k');
        plot([x(i)+xe(i) x(i)+xe(i)],[y(i)-xticklen y(i)+xticklen],'k');
      end;
    else
      if (xe(i,1)>0)
        plot([x(i)-xe(i,1) x(i)],[y(i) y(i)],'k');
        plot([x(i)-xe(i,1) x(i)-xe(i,1)],[y(i)-xticklen y(i)+xticklen],'k');
      end;
      if (xe(i,2)>0)
        plot([x(i)+xe(i,2) x(i)],[y(i) y(i)],'k');
        plot([x(i)+xe(i,2) x(i)+xe(i,2)],[y(i)-xticklen y(i)+xticklen],'k');
      end;
    end;

    if (isvector(ye))
      if (ye(i)>0)
        plot([x(i) x(i)],[y(i)-ye(i) y(i)+ye(i)],'k');
        plot([x(i)-yticklen x(i)+yticklen],[y(i)-ye(i) y(i)-ye(i)],'k');
        plot([x(i)-yticklen x(i)+yticklen],[y(i)+ye(i) y(i)+ye(i)],'k');
      end;
    else
      if (ye(i,1)>0)
        plot([x(i) x(i)],[y(i)-ye(i,1) y(i)],'k');
        plot([x(i)-yticklen x(i)+yticklen],[y(i)-ye(i,1) y(i)-ye(i,1)],'k');
      end;
      if (ye(i,2)>0)
        plot([x(i) x(i)],[y(i)+ye(i,2) y(i)],'k');
        plot([x(i)-yticklen x(i)+yticklen],[y(i)+ye(i,2) y(i)+ye(i,2)],'k');
      end;
    end;
  end;
  hold off;
  putbnd([x-xe(:,1); x+xe(:,2)],[y-ye(:,1); y+ye(:,2)]);

  return;
