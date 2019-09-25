% PUTREGRLINE: Add a predictive or major-axis regression line to the current 
%              plot.
%
%     Usage: b = putregrline(x,y,{kind},{w})
%
%         x,y =   abscissa and ordinate values of current plot.
%         kind =  optional flag indicating kind of regression:
%                   0 = predictive regression [default];
%                   1 = major-axis regression.
%         w =     optional vector of weights for weighted predictive regression
%                   [default = null].
%         ---------------------------------------------------------------------
%         b =     regression coefficients.
%

% RE Strauss, 4/10/01
%   11/9/01 - corrected error with implementation of 'kind' flag.
%   5/29/02 - change function name from 'addregrline'.
%   11/26/02 - return the regression coefficients.
%   7/24/03 - added option of weighted predictive regression.

function b = putregrline(x,y,kind,w)
  if (nargin < 3) kind = []; end;
  if (nargin < 4) w = []; end;

  if (isempty(kind))
    kind = 0;
  end;

  switch (kind)
    case 0,                               % Major-axis regression
      [b,stats,pred] = linregr(x,y,w);

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
