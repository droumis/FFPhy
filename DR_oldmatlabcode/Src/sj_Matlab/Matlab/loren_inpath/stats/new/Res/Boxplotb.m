% BOXPLOTB: Plots single bar for a box plot, assuming that the current plot is being 
%           held ('hold on');
%
%     Usage: [lm,rm] = boxplotb(y,x,width,{nparm})
%
%           y =     vector containing distribution of values for a single group
%           x =     abscissa value at which boxplot is to be plotted
%           width = width of boxplot
%           nparm = optional boolean flag indicating summary statistics to be plotted:
%                     0 = mean, stdev, range [default]
%                     1 = quartiles
%           --------------------------------------------------------------------------
%           lm =    [1 x 2] vector of coordinates of left edge of mean/median bar.
%           rm =    [1 x 2] vector of coordinates of right edge of mean/median bar.
%

% RE Strauss, 5/6/98
%   9/7/99 - changed plot colors for Matlab v5.

function [lm,rm] = boxplotb(y,x,width,nparm)
  if (nargin < 4) nparm = []; end;

  if (isempty(nparm))
    nparm = 0;
  end;

  ystat = zeros(5,1);                         % Box positions, bottom to top
  ystat(1) = min(y);
  ystat(5) = max(y);

  if (nparm)
    p = [25,50,75];
    ystat(2:4) = prctile(y,p)';
  else
    ym = mean(y);
    ys = std(y);
    ystat(2) = ym - ys;
    ystat(3) = ym;
    ystat(4) = ym + ys;
  end;

  xl = x - width/2;
  xh = xl + width;

  for i = 1:5
    plot([xl xh],[ystat(i) ystat(i)],'k');    % Crossbars
  end;
  plot([xl xl],[ystat(2) ystat(4)],'k');
  plot([xh xh],[ystat(2) ystat(4)],'k');
  plot([x x],[ystat(1) ystat(2)],'k');
  plot([x x],[ystat(4) ystat(5)],'k');

  lm = [xl ystat(3)];
  rm = [xh ystat(3)];

  return;
