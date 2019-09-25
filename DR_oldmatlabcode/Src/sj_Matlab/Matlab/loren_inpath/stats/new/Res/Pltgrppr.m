% PLTGRPPR: Plot regression lines for function plotgrps()
%

% RE Strauss, 5/31/99
%   9/2/99 - changed plot colors for Matlab v5.

function XY = pltgrppr(x,y,XY)
  xdev = x-mean(x);
  ydev = y-mean(y);
  slope = (xdev'*ydev) / (xdev'*xdev);
  intcpt = mean(y) - slope*mean(x);

  xmargin = 0.05*range(x);
  ymargin = 0.05*range(y);

  xmin = min(x);
  xmax = max(x);
  ymin = min(y);
  ymax = max(y);

  xlow =  xmin - xmargin;
  xhigh = xmax + xmargin;

  ylow = xlow * slope + intcpt;
  yhigh = xhigh * slope + intcpt;
        
  if (yhigh > ymax+ymargin)
          yhigh = ymax + ymargin;
          xhigh = (yhigh - intcpt) / slope;
  end;

  if (ylow < ymin-ymargin)
          ylow = ymin - ymargin;
          xlow = (ylow - intcpt) / slope;
  end;

  plot([xlow,xhigh],[ylow,yhigh],'k');

  XY = [XY; xlow ylow; xhigh yhigh];  % Extend axis ranges

  return;

