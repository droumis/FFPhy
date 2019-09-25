% ADDREGRLINE: Add a predictive or major-axis regression line to the current 
%              plot.
%
%     Usage: addregrline(x,y,kind)
%
%         x,y =   abscissa and ordinate values of current plot.
%         kind =  optional flag indicating kind of regression:
%                   0 = predictive regression [default]
%                   1 = major-axis regression
%

% RE Strauss, 4/10/01
%   11/9/01 - corrected error with implementation of 'kind' flag.

function addregrline(x,y,kind)
  if (nargin < 3) kind = []; end;

  if (isempty(kind))
    kind = 0;
  end;

  switch (kind)
    case 0,                               % Major-axis regression
      [b,stats,pred] = linregr(x,y);

    case 1,                               % Predictive regression
      b = majaxis(x,y);
      pred = [x ones(size(x))]*b';
  end;

  [xmin,imin] = min(x);
  [xmax,imax] = max(x);
  hold on;
  plot([x([imin imax])],[pred([imin imax])],'k');
  hold off;

  return;
