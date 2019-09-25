% NORMDRAW: Plots one or more normal distributions
%
%     Usage: normdraw(m,s,{bound})
%
%           m =     vector of means.
%           s =     corresponding vector of standard deviations.
%           bound = 2-element vector of lower and upper bounds of plot
%                     [default = min(m)-min(3s), max(m)+max(3s)].
%

% RE Strauss, 3/5/99
%   9/7/99 -  changed plot colors for Matlab v5.
%   2/24/00 - added optional bounds; remove creation of new figure window.
%   5/5/01 -  changed default bounds.

function normdraw(m,s,bound)
  if (nargin < 3) bound = []; end;

  if (isempty(bound))
    bound = [min(m)-max(3*s) max(m)+max(3*s)];
  end;
  if (bound(2)<=bound(1))
    error('  NORMDRAW: invalid bounds');
  end;

  x = linspace(bound(1),bound(2));
  maxfx = 0;

  fx = normpdf(x,m(1),s(1));
  if (max(fx)>maxfx)
    maxfx = max(fx);
  end;

  plot(x,fx,'k');
  hold on;

  for i = 2:length(m)
    fx = normpdf(x,m(i),s(i));
    if (max(fx)>maxfx)
      maxfx = max(fx);
    end;
    plot(x,fx,'k');
  end;

  hold off;
  delta = 0.05*(bound(2)-bound(1));
  axis([bound(1)-delta, bound(2)+delta, 0, 1.05*maxfx]);
  putxlab('X');
  putylab('f(X)');

  return;