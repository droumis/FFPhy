% BILINEARF:  Objective function for bilinear().  Fits reparameterized model of 
%             Chappell (1989).
%
%     Usage: [mse,b,crd] = bilinearf(t,x,y)
%
%         t =   scalar indicating the abscissa of the change-point.
%         x =   vector of independent-variable values. 
%         y =   corresponding vector of dependent-variable values.
%         ---------------------------------------------------------------------
%         mse = mean residual sum-of-squares.
%         b =   reparameterized parameters [a1,a2,b1,b2], where a_i are the 
%                 intercepts of the two linear segments and b_i are the slopes.
%         crd = [3 x 2] matrix of coordinates of endpoints of line segments.
%

% Chappell, R. 1989. Fitting bent lines to data, with applications to 
%   allometry.  J. Theor. Biol. 138:235-256.
%   [NB error on p. 239: min(0,x-t) should be max(0,x-t)].

% RE Strauss, 10/8/00
  
function [mse,b,crd] = bilinearf(t,x,y)
  get_crd = 0;
  if (nargout > 2)
    get_crd = 1;
  end;

  x = x(:);
  y = y(:);
  n = length(x);

  if (t<=min(x) | t>=max(x))
    error('  BILINEARF: t out of range of data.');
  end;

  tt = t*ones(n,1);

  u = min(x,tt);
  v = max(zeros(n,1),x-tt);

  X = [u v ones(n,1)];                    % Multiple regression
  b = inv(X'*X) * X'*y;

  a1 = b(3);
  a2 = b(3)+(b(1)-b(2))*t;
  b1 = b(1);
  b2 = b(2);

  b = [a1 a2 b1 b2];

  e = y - (a1 + b1*u + b2*v);
  mse = (e'*e)/(n-4);

  if (get_crd)
    u = [min(x) t t];
    v = [0 0 max(x)-t];
    crd = [min(x) t max(x); (a1 + b1*u + b2*v)]';
  end;

  return;
