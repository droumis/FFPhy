% CIRCCRDS: Returns plotting coordinatess for a circle specified by its circle 
%           and radius.
%
%     Usage: crds = circcrds({radius},{center},{npts})
%
%           radius = radius of circle [default = 1].
%           center = [1 x 2] row vector of point coordinates of fitted center
%                     [default = origin].
%           npts =   optional number of plotting points to be returned 
%                     [default = 100].
%           -----------------------------------------------------------------
%           crds =   [npts x 2] matrix of point coordinates.
%

% RE Strauss, 2/24/98
%   10/7/99 - changed sequence of input arguments.
%   3/19/00 - fixed error to return column vectors.
%   6/29/01 - added error messages.

function crds = circcrds(radius,center,npts)
  if (nargin < 1) radius = []; end;
  if (nargin < 2) center = []; end;
  if (nargin < 3) npts = []; end;

  if (isempty(radius))
    radius = 1;
  end;
  if (isempty(center))
    center = [0 0];
  end;
  if (isempty(npts))
    npts = 100;
  end;

  center = center(:);

  if (~isscalar(radius))
    error('  CIRCCRDS: radius must be a scalar');
  end;
  if (length(center)~=2)
    error('  CIRCCRDS: center must be vector of length 2');
  end;
  if (~isscalar(npts))
    error('  CIRCCRDS: npts must be a scalar');
  end;

  theta = linspace(0,2*pi,npts)';
  x = radius .* cos(theta) + center(1);
  y = radius .* sin(theta) + center(2);

  crds = [x y];

  return;


