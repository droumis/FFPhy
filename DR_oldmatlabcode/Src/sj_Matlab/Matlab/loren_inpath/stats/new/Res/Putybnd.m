% PUTYBND:  Changes the [min,max] settings for the y-axis, leaving the 
%           y-axis settings unchanged.
%
%     Syntax: putybnd(ymin,ymax)
%

% RE Strauss, 12/8/96
%   2/2/00 - changed handling of input arguments.

function putybnd(ymin,ymax)
  if (nargin < 2) ymax = []; end;

  if (isempty(ymax) & max(size(ymin))==2)
    ymax = ymin(2);
    ymin = ymin(1);
  end;

  if (~isscalar(ymin) | ~isscalar(ymax))
    error('  PUTXBND: invalid axis settings');
  end;

  v = axis;
  v(3) = ymin;
  v(4) = ymax;
  axis(v);

  return;
