% BILINTRANSF:  Objective function for bilintrans().  Fits Koops & Grossman 
%               (1993) model.
%
%     Usage: [mse,crds] = bilintransf(p,x,y)
%
%         b =   [1 x 5] vector of parameters: [a,b1,b2,c,r].
%         x =   vector of independent-variable values. 
%         y =   corresponding vector of dependent-variable values.
%         ---------------------------------------------------------------------
%         mse = mean residual sum-of-squares.
%         crds = [3 x 2] matrix of coordinates of endpoints of line segments.
%

% Koops, W.J. & M. Grossman. 1993. Multiphasic allometry. Growth, Development & 
%   Aging 57:183-192.

% RE Strauss, 10/9/00

function [mse,crds] = bilintransf(p,x,y)
  get_crds = 0;
  if (nargout > 1)
    get_crds = 1;
  end;

  x = x(:);
  y = y(:);
  n = length(x);

  [a,b1,b2,c,r] = extrcols(p);          % Extract parameters from vector

  if (c<=min(x) | c>=max(x))
%    error('  BILINTRANSF: c out of range of data.');
    mse = 1e6;
    crds = [];
  end;

  ypred = a + b1*x - (b1-b2)*r*log(1+exp((x-c)/r));
  e = y - ypred;

  mse = (e'*e)/(n-5);

  if (get_crds)
    cx = linspace(min(x),max(x))';
    cy = a + b1*cx - (b1-b2)*r*log(1+exp((cx-c)/r));
    crds = [cx cy];
  end;

  return;

