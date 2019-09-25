% PEARSONN: Pearson transitional curve: normal
%
%     Usage: [xx,yy,mse] = pearsonn(x,fx,mu,xmin,xmax)
%
%         x =       midpoints of histogram bars.
%         fx =      heights of histogram bars.
%         mu =      vector of moments.
%         xmin =    lower extent of function on axis.
%         xmax =    upper extent of function on axis.
%         ----------------------------------------------------------
%         x,y =     coordinates of fitted function.
%         mse =     mean squared deviation between fx and predicted.
%

% RE Strauss, 3/4/01

function [xx,yy,mse] = pearsonn(x,fx,mu,xmin,xmax)
  y0 = 1;
  yy = y0*normpdf(x,mu(1),sqrt(mu(2)));
  y0 = y0/sum(yy);
  yy = y0*normpdf(x,mu(1),sqrt(mu(2)));
  mse = abs(mean((yy-fx).^2));

  xx = linspace(xmin,xmax)';
  yy = y0*normpdf(xx,mu(1),sqrt(mu(2)));

  return;
