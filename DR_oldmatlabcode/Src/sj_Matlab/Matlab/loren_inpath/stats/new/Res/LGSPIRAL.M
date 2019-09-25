% LGSPIRAL: Generates figure of a logarithmic spiral, based on the two 
%           parameters alpha and b.  Alpha should be in the approximate 
%           range of 1.2-1.6 radians.  If alpha = pi/2 = 1.57 radians, 
%           the spiral degerates into a circle; if alpha is greater than 
%           pi/2, the spiral spirals inwards rather than outwards.  
%           b is a scaling factor. 
%
%     Usage: lgspiral(alpha,b,whorls)
%
%           alpha =  angular displacement of spiral tangent from radial 
%                      line (shape) [default = 1.4 radians].
%           b =      scale factor (size) [default = 1].
%           whorls = number of whorls to be plotted [default = 3].
%

% RE Strauss, 9/14/97
%   9/7/99 - changed plot colors for Matlab v5.

function lgspiral(alpha,b,whorls)
  if (nargin < 1) alpha = []; end;
  if (nargin < 2) b = []; end;
  if (nargin < 3) whorls = []; end;

  if (isempty(alpha))
    alpha = 1.4;
  end;
  if (isempty(b))
    b = 1;
  end;
  if (isempty(whorls))
    whorls = 3;
  end;

  a = cot(alpha);
  theta = linspace(0,2*pi*whorls,100*whorls);
  r = b*exp(a*theta);
  [x,y] = polarcrd(r,theta,1);

  figure;
  plot(x,y,'k-');
  v = putbnd(x,y);
%  axis('off');
  axis('equal');

  return;  
