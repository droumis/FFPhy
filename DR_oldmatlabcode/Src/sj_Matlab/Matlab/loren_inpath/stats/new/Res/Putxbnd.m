% PUTXBND:  Changes the [min,max] settings for the x-axis, leaving the 
%           y-axis settings unchanged.
%
%     Syntax: putxbnd(xmin,xmax)
%

% RE Strauss, 12/8/96
%   2/2/00 - changed handling of input arguments.

function putxbnd(xmin,xmax)
  if (nargin < 2) xmax = []; end;

  if (isempty(xmax) & max(size(xmin))==2)
    xmax = xmin(2);
    xmin = xmin(1);
  end;

  if (~isscalar(xmin) | ~isscalar(xmax))
    error('  PUTXBND: invalid axis settings');
  end;

  v = axis;
  v(1) = xmin;
  v(2) = xmax;
  axis(v);

  return;
