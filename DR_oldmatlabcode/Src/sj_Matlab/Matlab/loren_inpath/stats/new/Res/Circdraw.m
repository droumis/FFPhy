% CIRCDRAW: Plots a circle for a given radius, r.  Assumes that other plot 
%             characteristics (e.g., square plot) have been set by calling 
%             function.
%
%   Usage: circdraw(r)
%
%     r = optional radius [default = 1].
%

% RE Strauss, 3/25/98
%   8/21/99 - set default radius; changed plot colors for Matlab v5.

function circdraw(r)
  if (nargin < 1)
    r = 1;
  end;

  t = pi-(-1:.01:1)'*pi;                % Draw circle
  x = [r*cos(t)];
  y = [r*sin(t)];
  plot(x,y,'k');

  return;

