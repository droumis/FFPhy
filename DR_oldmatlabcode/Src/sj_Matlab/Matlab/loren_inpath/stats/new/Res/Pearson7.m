% PEARSON7: Pearson Type VII curve
%
%     Usage: [xx,yy,mse] = pearson7(x,fx,mu,b1,b2,xmin,xmax)
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

function [xx,yy,mse] = pearson7(x,fx,mu,b1,b2,xmin,xmax)
  x = x - mu(1);

  m = (5*b2-9)/(2*b2-6);
  a2 = (2*mu(2)*b2)/(b2-3);

  y0 = 1;
  yy = y0*(1-x.^2./a2).^(-m);
  y0 = y0/sum(yy);
  yy = y0*(1-x.^2./a2).^(-m);
  mse = abs(mean((yy-fx).^2));

  xx = linspace(xmin-mu(1),xmax-mu(1))';
  yy = y0*(1-xx.^2./a2).^(-m);

  xx = xx+mu(1);

  return;
