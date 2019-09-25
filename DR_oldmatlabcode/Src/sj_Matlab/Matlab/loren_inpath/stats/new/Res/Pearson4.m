% PEARSON4: Pearson Type IV curve
%
%     Usage: [xx,yy,mse] = pearson4(x,fx,mu,b1,b2,xmin,xmax)
%
%         x =       midpoints of histogram bars.
%         fx =      heights of histogram bars.
%         mu =      vector of moments.
%         b1,b2 =   beta values.
%         xmin =    lower extent of function on axis.
%         xmax =    upper extent of function on axis.
%         ----------------------------------------------------------
%         x,y =     coordinates of fitted function.
%         mse =     mean squared deviation between fx and predicted.
%

% RE Strauss, 3/4/01

function [xx,yy,mse] = pearson4(x,fx,mu,b1,b2,xmin,xmax)
  x = x - mu(1);

  r = 6*(b2-b1-1)/(2*b2-3*b1-6);
  rad = sqrt(16*(r-1)-b1*(r-2)^2);
  m = 0.5*(r+2);
  a = sqrt(mu(2)/16)*rad;
  v = abs((-r*(r-2)*sqrt(b1))/rad) * -sign(mu(3));

  y0 = 1;
  yy = y0 .* (1+(x./a-v/r).^2).^(-m) .* exp(-v*atan(x./a-v/r));
  y0 = y0/sum(yy);

  yy = y0 .* (1+(x./a-v/r).^2).^(-m) .* exp(-v*atan(x./a-v/r));
  mse = abs(mean((yy-fx).^2));

  xx = linspace(xmin-mu(1),xmax-mu(1))';
  yy = y0 .* (1+(xx./a-v/r).^2).^(-m) .* exp(-v*atan(xx./a-v/r));
  xx = xx + mu(1);

  return;
