% PLOTBERN: Plots a relative Bernoulli distribution (X={0,1}), given absolute or 
%           relative frequencies of 0's and 1's.
%
%     Usage: plotbern(f0,f1,{bwidth})
%
%           f0 =     f(0)
%           f1 =     f(1)
%           bwidth = optional bar width, in scale units [default = 0.3].
%

% RE Strauss, 6/7/97
%   11/25/99 - changed plot colors for Matlab v5 and scaling of abscissa.

function plotbern(f0,f1,bwidth)
  if (nargin < 3) bwidth = []; end;

  if (isempty(bwidth))
    bwidth = 0.3;
  end;

  s = f0+f1;
  if (s>1)
    f0 = f0/s;
    f1 = f1/s;
  end;

  w = bwidth;      % Bar width

  x = [-w/2 -w/2 w/2 w/2];
  y = [0 f0 f0 0];
  plot(x,y,'k');
  hold on;
  fill(x,y,'k');

  x = [1-w/2 1-w/2 1+w/2 1+w/2];
  y = [0 f1 f1 0];
  plot(x,y,'k');
  fill(x,y,'k');

  v = [-w 1+w 0 1.05];
  axis(v);
  puttick([0 1]);
  putylab('Relative frequency');
  hold off;

  return;
